use ark_bls12_381::G1Affine;
use ark_ec::ProjectiveCurve;
use ark_std::Zero;

use crate::{
    batch_adder::BatchAdder,
    bitmap::Bitmap,
    glv::endomorphism,
    types::{G1Projective, GROUP_SIZE, GROUP_SIZE_LOG2},
};

pub struct BucketMSM {
    num_windows: u32,
    window_bits: u32,
    bucket_bits: u32,
    max_batch_cnt: u32,          // max slices allowed in a batch
    max_collision_cnt: u32,
    buckets: Vec<G1Affine>, // size (num_windows << window_bits) * 2

    // current batch state
    bitmap: Bitmap,
    batch_buckets_and_points: Vec<(u32, u32)>,
    collision_buckets_and_points: Vec<(u32, G1Affine)>,
    cur_points: Vec<G1Affine>, // points of current batch, size batch_size

    // batch affine adder
    batch_adder: BatchAdder,
}

impl BucketMSM {
    pub fn new(
        scalar_bits: u32,
        window_bits: u32,
        max_batch_cnt: u32,     // default: 4096
        max_collision_cnt: u32, // default: 128
    ) -> BucketMSM {
        let num_windows = (scalar_bits + window_bits - 1) / window_bits;
        let batch_size = std::cmp::max(8192, max_batch_cnt);
        let bucket_bits = window_bits - 1; // half buckets needed because of signed-bucket-index
        let bucket_size = num_windows << bucket_bits;
        // size of batch_adder will be the max of batch_size and num_windows * groups per window
        let batch_adder_size = std::cmp::max(batch_size, bucket_size >> GROUP_SIZE_LOG2);

        BucketMSM {
            num_windows,
            window_bits,
            bucket_bits,
            max_batch_cnt,
            max_collision_cnt,
            buckets: vec![G1Affine::zero(); bucket_size as usize],

            bitmap: Bitmap::new(bucket_size as usize / 32),
            batch_buckets_and_points: Vec::with_capacity(batch_size as usize),
            collision_buckets_and_points: Vec::with_capacity(max_collision_cnt as usize),
            cur_points: vec![G1Affine::zero(); batch_size  as usize],

            batch_adder: BatchAdder::new(batch_adder_size as usize),
        }
    }

    pub fn process_point_and_slices_glv(
            &mut self, point: &G1Affine,
            normal_slices: &Vec<u32>,
            phi_slices: &Vec<u32>,
            is_neg_scalar: bool,
            is_neg_normal: bool) {
        assert!(
            self.num_windows as usize == normal_slices.len() && normal_slices.len() == phi_slices.len(),
            "slice len check failed: normal_slices {}, phi_slices {}, num_windows {}",
            normal_slices.len(), phi_slices.len(), self.num_windows
        );

        let mut p = *point; // copy

        if is_neg_scalar {p.y = -p.y};
        if is_neg_normal {p.y = -p.y};

        self.cur_points.push(p.clone());
        for win in 0..normal_slices.len() {
            if (normal_slices[win] as i32) > 0 {
                let bucket_id = (win << self.bucket_bits) as u32 + normal_slices[win] - 1;
                self._process_slices(bucket_id, self.cur_points.len() as u32 - 1);
            }
        }

        p.y = -p.y;

        self.cur_points.push(p.clone());
        for win in 0..normal_slices.len() {
            if (normal_slices[win] as i32) < 0 {
                let slice = normal_slices[win] & 0x7FFFFFFF;
                if slice > 0 {
                    let bucket_id = (win << self.bucket_bits) as u32 + slice - 1;
                    self._process_slices(bucket_id, self.cur_points.len() as u32 - 1);
                }
            }
        }

        // process phi slices
        p.y = -p.y;
        if is_neg_normal {p.y = -p.y;}
        endomorphism(&mut p);

        self.cur_points.push(p.clone());
        for win in 0..phi_slices.len() {
            if (phi_slices[win] as i32) > 0 {
                let bucket_id = (win << self.bucket_bits) as u32 + phi_slices[win] - 1;
                self._process_slices(bucket_id, self.cur_points.len() as u32 - 1);
            }
        }

        p.y = -p.y;

        self.cur_points.push(p.clone());
        for win in 0..phi_slices.len() {
            if (phi_slices[win] as i32) < 0 {
                let slice = phi_slices[win] & 0x7FFFFFFF;
                if slice > 0 {
                    let bucket_id = (win << self.bucket_bits) as u32 + slice - 1;
                    self._process_slices(bucket_id, self.cur_points.len() as u32 - 1);
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn process_point_and_slices(&mut self, point: &G1Affine, slices: &Vec<u32>) {
        assert!(
            self.num_windows as usize == slices.len(),
            "slices.len() {} should equal num_windows {}",
            slices.len(), self.num_windows
        );

        self.cur_points.push(point.clone());
        for win in 0..slices.len() {
            if (slices[win] as i32) > 0 {
                let bucket_id = (win << self.bucket_bits) as u32 + slices[win] - 1; // skip slice == 0
                self._process_slices(bucket_id, self.cur_points.len() as u32 - 1);
            }
        }

        let mut neg_p = *point;
        neg_p.y = -neg_p.y;

        self.cur_points.push(neg_p);
        for win in 0..slices.len() {
            if (slices[win] as i32) < 0 {
                let slice = slices[win] & 0x7FFFFFFF;
                if slice > 0 {
                    let bucket_id = (win << self.bucket_bits) as u32 + slice - 1; // skip slice == 0
                    self._process_slices(bucket_id, self.cur_points.len() as u32 - 1);
                }
            }
        }
    }

    pub fn process_complete(&mut self) {
        self._process_batch();
        while !(self.collision_buckets_and_points.is_empty() && self.batch_buckets_and_points.is_empty()) {
            self._process_batch();
        }
    }

    fn _process_slices(&mut self, bucket_id: u32, point_idx: u32) {
        if !self.bitmap.test_and_set(bucket_id) {
            // if no collision found, add point to current batch
            self.batch_buckets_and_points.push((bucket_id, point_idx));
        } else {
            self.collision_buckets_and_points.push((bucket_id, self.cur_points[point_idx as usize]));
        }

        if self.collision_buckets_and_points.len() as u32 >= self.max_collision_cnt ||
                self.batch_buckets_and_points.len() as u32 >= self.max_batch_cnt {
            self._process_batch();
        }
    }

    fn _process_batch(&mut self) {
        if self.batch_buckets_and_points.is_empty() {
            return;
        }
        // batch addition
        let (bucket_ids, point_idxs): (Vec<u32>, Vec<u32>) = self.batch_buckets_and_points
            .iter().map(|(b, p)| (*b, *p)).unzip();
        self.batch_adder.batch_add_indexed(&mut self.buckets, &bucket_ids, &self.cur_points, &point_idxs);
        // clean up current batch
        self.bitmap.clear();
        self.batch_buckets_and_points.clear();
        // memorize the last point which is the current processing point and we need to
        // push it back to the cur_points list since we're processing slices in a for loop
        let slicing_point = self.cur_points.pop();
        self.cur_points.clear();

        let mut next_pos = 0;
        for i in 0..self.collision_buckets_and_points.len() {
            let (bucket_id, point) = self.collision_buckets_and_points[i];
            if self.bitmap.test_and_set(bucket_id) {
                // collision found
                self.collision_buckets_and_points.swap(next_pos, i);
                next_pos += 1;
            } else {
                self.batch_buckets_and_points.push((bucket_id, self.cur_points.len() as u32));
                self.cur_points.push(point);
            }
        }
        self.collision_buckets_and_points.truncate(next_pos);
        self.cur_points.push(slicing_point.unwrap());
    }

    pub fn batch_reduce(&mut self) -> G1Projective {
        let window_starts: Vec<_> = (0..self.num_windows as usize).collect();
        let num_groups = (self.num_windows as usize) << (self.bucket_bits as usize - GROUP_SIZE_LOG2);
        let mut running_sums: Vec<_> = vec![G1Affine::zero(); num_groups];
        let mut sum_of_sums: Vec<_> = vec![G1Affine::zero(); num_groups];

        // calculate running sum and sum of sum for each group
        for i in (0..GROUP_SIZE).rev() {
            // running sum
            self.batch_adder.batch_add_step_n(&mut running_sums,
                                              1,
                                              &self.buckets[i..],
                                              GROUP_SIZE,
                                              num_groups);
            // sum of sum
            self.batch_adder.batch_add(&mut sum_of_sums, &running_sums);
        }

        let sum_by_window: Vec<G1Projective> = ark_std::cfg_into_iter!(window_starts)
            .map(|w_start| {
                let group_start = w_start << (self.bucket_bits as usize - GROUP_SIZE_LOG2);
                let group_end = (w_start + 1) << (self.bucket_bits as usize - GROUP_SIZE_LOG2);
                self.inner_window_reduce(&running_sums[group_start..group_end],
                                         &sum_of_sums[group_start..group_end])
            }).collect();

        return self.intra_window_reduce(&sum_by_window);
    }

    fn inner_window_reduce(&mut self, running_sums: &[G1Affine], sum_of_sums: &[G1Affine]) -> G1Projective {
        return self.calc_sum_of_sum_total(sum_of_sums) + self.calc_running_sum_total(running_sums);
    }

    fn calc_running_sum_total(&mut self, running_sums: &[G1Affine]) -> G1Projective {
        let mut running_sum_total = G1Projective::zero();
        for i in 1..running_sums.len() {
            for _ in 0..i {
                running_sum_total.add_assign_mixed(&running_sums[i]);
            }
        }

        for _ in 0..GROUP_SIZE_LOG2 {
            running_sum_total.double_in_place();
        }
        return running_sum_total;
    }

    fn calc_sum_of_sum_total(&mut self, sum_of_sums: &[G1Affine]) -> G1Projective {
        let mut sum = G1Projective::zero();
        sum_of_sums.iter().for_each(|p| sum.add_assign_mixed(p));
        return sum;
    }

    fn intra_window_reduce(&mut self, window_sums: &Vec<G1Projective>) -> G1Projective {
        // We store the sum for the lowest window.
        let lowest = *window_sums.first().unwrap();

        // We're traversing windows from high to low.
        lowest
        + &window_sums[1..]
            .iter()
            .rev()
            .fold(G1Projective::zero(), |mut total, sum_i| {
                total += sum_i;
                for _ in 0..self.window_bits {
                    total.double_in_place();
                }
                total
            })
    }
}

#[cfg(test)]
mod bucket_msm_tests {
    use super::*;
    use ark_bls12_381::G1Affine;
    use ark_std::UniformRand;
    use ark_ec::AffineCurve;

    #[test]
    fn test_process_point_and_slices_deal_two_points() {
        let window_bits = 15u32;
        let mut bucket_msm = BucketMSM::new(30u32, window_bits, 128u32, 4096u32);
        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let q_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p = G1Affine::from(p_prj);
        let q = G1Affine::from(q_prj);

        bucket_msm.process_point_and_slices(&p, &vec![1u32, 3u32]);
        bucket_msm.process_point_and_slices(&q, &vec![2u32, 3u32]);
        bucket_msm.process_complete();
        assert_eq!(bucket_msm.buckets[0], p);
        assert_eq!(bucket_msm.buckets[1], q);
        assert_eq!(bucket_msm.buckets[2 + (1 << bucket_msm.bucket_bits)], p + q);
    }

    #[test]
    fn test_process_point_and_slices_deal_three_points() {
        let window_bits = 15u32;
        let mut bucket_msm = BucketMSM::new(45u32, window_bits, 128u32, 4096u32);
        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let q_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let r_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p = G1Affine::from(p_prj);
        let q = G1Affine::from(q_prj);
        let r = G1Affine::from(r_prj);

        bucket_msm.process_point_and_slices(&p, &vec![1u32, 3u32, 4u32]);
        bucket_msm.process_point_and_slices(&q, &vec![2u32, 3u32, 4u32]);
        bucket_msm.process_point_and_slices(&r, &vec![2u32, 3u32, 5u32]);
        bucket_msm.process_complete();
        assert_eq!(bucket_msm.buckets[0], p);
        assert_eq!(bucket_msm.buckets[1], q + r);
        assert_eq!(
            bucket_msm.buckets[2 + (1 << bucket_msm.bucket_bits)],
            p + q + r
        );
        assert_eq!(bucket_msm.buckets[3 + (2 << bucket_msm.bucket_bits)], p + q);
        assert_eq!(bucket_msm.buckets[4 + (2 << bucket_msm.bucket_bits)], r);
    }

    #[test]
    fn test_process_point_and_slices_glv_deal_two_points() {
        let window_bits = 15u32;
        let mut bucket_msm = BucketMSM::new(30u32, window_bits, 128u32, 4096u32);
        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let q_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let mut p = G1Affine::from(p_prj);
        let mut q = G1Affine::from(q_prj);

        bucket_msm.process_point_and_slices_glv(&p, &vec![1u32, 3u32], &vec![4u32, 6u32], false, false);
        bucket_msm.process_point_and_slices_glv(&q, &vec![2u32, 3u32], &vec![5u32, 6u32], false, false);
        bucket_msm.process_complete();
        assert_eq!(bucket_msm.buckets[0], p);
        assert_eq!(bucket_msm.buckets[1], q);
        assert_eq!(bucket_msm.buckets[2 + (1 << bucket_msm.bucket_bits)], p + q);

        endomorphism(&mut p);
        endomorphism(&mut q);
        assert_eq!(bucket_msm.buckets[3], p);
        assert_eq!(bucket_msm.buckets[4], q);
        assert_eq!(bucket_msm.buckets[5 + (1 << bucket_msm.bucket_bits)], p + q);
    }
}
