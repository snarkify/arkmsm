use ark_bls12_381::G1Affine;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_std::{Zero};

use crate::{bitmap::Bitmap, list, types::G1Projective, batch_adder::BatchAdder};

const GROUP_SIZE_LOG2: usize = 6;
const GROUP_SIZE: usize = 1 << GROUP_SIZE_LOG2;

pub struct BucketMSM<G: AffineCurve> {
    num_windows: u32,
    window_bits: u32,
    bucket_bits: u32,

    buckets: Vec<G1Affine>, // size (num_windows << window_bits) * 2
    results: G,             // output

    // current batch state
    cur_points: Vec<G1Affine>, // points of current batch, size batch_size + 2 x max_collision_cnt
    cur_points_cnt: u32,
    cur_bkt_pnt_pairs: Vec<u64>, // encoded as <point_idx, bucket> pair, size max_batch_cnt
    cur_batch_cnt: u32,          // number of slices to be processed in current batch
    max_batch_cnt: u32,          // max slices allowed in a batch

    // batch affine adder
    batch_adder: BatchAdder,

    // collision and pending slices state
    //
    // a collision occurs when adding a point to bucket and found is already a
    // point in the same batch has been assigned to the same bucket
    bitmap: Bitmap,
    collision_cnt: u32,
    max_collision_cnt: u32, // max collisions allowed per batch

    // worst-case: all slices, except for the first slice, from the
    // unprocessed_list were added to the processing_list, and a new set of
    // collisions occur, causing the unprocessed_list to be filled again. To
    // accommodate both lists, we must allocate 2x max_collision_cnt
    collision_list_nodes: Vec<u64>, // u32 encoded as <point index to cur_points, next> pair, size 2 x max_collision_cnt
    available_list: u32,            // free entries in collision_list_nodes
    processing_list: u32,           // entries to be processed in the next batch
    unprocessed_list: u32,          // pending processing entries
}

impl<G: AffineCurve> BucketMSM<G> {
    pub fn new(
        scalar_bits: u32,
        window_bits: u32,
        max_batch_cnt: u32,     // default: 4096
        max_collision_cnt: u32, // default: 128
    ) -> BucketMSM<G> {
        let num_windows = (scalar_bits + window_bits - 1) / window_bits;
        let batch_size = std::cmp::max(8192, max_batch_cnt);
        let bucket_bits = window_bits - 1; // half buckets needed because of signed-bucket-index
        let bucket_size = num_windows << bucket_bits;
        // size of batch_adder will be the max of batch_size and num_windows * groups per window
        let batch_adder_size = std::cmp::max(batch_size, bucket_size >> GROUP_SIZE_LOG2);

        // link all the collision_list_nodes into a list
        let max_size = max_collision_cnt * 2;
        let mut _collision_list_nodes = vec![0u64; max_size as usize];
        for i in 0..(max_size - 1) {
            _collision_list_nodes[i as usize] = i as u64 + 1;
        }
        _collision_list_nodes[(max_size - 1) as usize] = list::NIL as u64;

        BucketMSM {
            num_windows,
            window_bits,
            bucket_bits,

            buckets: vec![G1Affine::zero(); bucket_size as usize],
            results: G::zero(),

            // 2 * max_collision_cnt for unprocessed and processing
            cur_points: vec![G1Affine::zero(); (batch_size + 2 * max_collision_cnt) as usize],
            cur_points_cnt: 0,
            cur_bkt_pnt_pairs: vec![0; batch_size as usize],
            cur_batch_cnt: 0,
            max_batch_cnt,

            batch_adder: BatchAdder::new(batch_adder_size as usize),

            bitmap: Bitmap::new(bucket_size as usize / 32),
            collision_cnt: 0,
            max_collision_cnt,

            collision_list_nodes: _collision_list_nodes,
            processing_list: list::make_list(list::NIL, list::NIL),
            unprocessed_list: list::make_list(list::NIL, list::NIL),
            available_list: list::make_list(0, max_size - 1),
        }
    }

    pub fn process_point_and_slices(&mut self, point: &G1Affine, slices: &Vec<u32>) {
        assert!(
            self.num_windows as usize == slices.len(),
            "slices.len() {} should equal num_windows {}",
            slices.len(),
            self.num_windows
        );

        // print!("point count: {}, num_windows {}, slices.len() {}\n", self.cur_points_cnt, self.num_windows, slices.len());
        // print!("slices {:?}\n", slices);

        self.cur_points[self.cur_points_cnt as usize] = *point; // copy
        self.cur_points_cnt += 1;
        for win in 0..slices.len() {
            if (slices[win] as i32) > 0 {
                let bucket_id = (win << self.bucket_bits) as u32 + slices[win] - 1; // skip slice == 0

                // print!("win: {}, win shift: {}, slice: {}, bucket_id {}\n", win, (win << self.bucket_bits), slices[win], bucket_id);
                self._process_slices(bucket_id, point);
            }
        }

        let mut neg_p = *point;
        neg_p.y = -neg_p.y;

        self.cur_points[self.cur_points_cnt as usize] = neg_p; // copy
        self.cur_points_cnt += 1;
        for win in 0..slices.len() {
            if (slices[win] as i32) < 0 {
                let slice = slices[win] & 0x7FFFFFFF;
                if slice > 0 {
                    let bucket_id = (win << self.bucket_bits) as u32 + slice - 1; // skip slice == 0

                    // print!("win: {}, win shift: {}, slice: {}, bucket_id {}\n", win, (win << self.bucket_bits), slices[win], bucket_id);
                    self._process_slices(bucket_id, &neg_p);
                }
            }
        }
    }

    pub fn process_complete(&mut self) {
        self._process_batch();
        while !list::is_empty(self.unprocessed_list) || !list::is_empty(self.processing_list) {
            self._process_batch();
        }
    }

    fn _process_slices(&mut self, bucket_id: u32, point: &G1Affine) {
        if !self.bitmap.test_and_set(bucket_id) {
            // if no collision found, add point to current batch

            let acc = &self.buckets[bucket_id as usize];
            // print!("bucket_id: {}, pair index: {}, point index {}\n", bucket_id, self.cur_batch_cnt, self.cur_points_cnt - 1);
            self.batch_adder.batch_add_phase_one(acc, point, self.cur_batch_cnt as usize);

            self.cur_bkt_pnt_pairs[self.cur_batch_cnt as usize] =
                ((bucket_id as u64) << 14) + self.cur_points_cnt as u64 - 1;
            self.cur_batch_cnt += 1;
        } else {
            // if collision, add point to unprocessed_list
            // free_slot index both collision_list_nodes and cur_points+max_collision_cnt
            //
            // processed_list:      ------------------------------->|
            // collision_list_nodes:            [<bucket_id, next>, <bucket_id, next>, ...]
            // cur_points: [ max_batch_cnt area |  point1,           point, ...           ]
            let free_node_index = list::pop(&self.collision_list_nodes, &mut self.available_list);
            // print!("bucket_id {}, free_node_index {}, collision_cnt {}, max_collision_cnt {}\n", bucket_id, free_node_index, self.collision_cnt, self.max_collision_cnt);
            assert!(free_node_index < 2 * self.max_collision_cnt);
            self.collision_list_nodes[free_node_index as usize] =
                ((bucket_id as u64) << 14) + 0x3FFF;
            list::enqueue(
                &mut self.collision_list_nodes,
                &mut self.unprocessed_list,
                free_node_index,
            );
            self.collision_cnt += 1;

            // cur_points: [max_batch_cnt points | 2*max_collisions points]
            self.cur_points[(self.max_batch_cnt + free_node_index) as usize] = *point; // copy
            // print!("Collision: bucket_id {}, point_idx {}\n", bucket_id, free_node_index);
        }

        if self.collision_cnt >= self.max_collision_cnt || self.cur_batch_cnt >= self.max_batch_cnt
        {
            self._process_batch();
        }
    }

    fn _process_batch(&mut self) {
        // print!("_process_batch: collision_cnt {} cur_batch_cnt {}\n", self.collision_cnt, self.cur_batch_cnt);

        // batch add phase two
        self.batch_adder.inverse();
        for i in (0..self.cur_batch_cnt).rev() {
            let bucket_point = self.cur_bkt_pnt_pairs[i as usize];
            let point = &self.cur_points[(bucket_point & 0x3FFF) as usize];
            let acc = &mut self.buckets[(bucket_point >> 14) as usize];
            // print!("phase_two: bucket_id: {}, pair index: {}, point_idx {}\n", bucket_point >> 14, i, bucket_point & 0x3FFF);
            // print!{"phase_two acc: {}\n", acc};
            // print!{"phase_two point: {}\n", point};
            self.batch_adder.batch_add_phase_two(acc, point, i as usize);
            // print!{"phase_two result: {}\n", acc};
        }

        // process collision
        self.bitmap.clear();

        list::append_all(
            &mut self.collision_list_nodes,
            &mut self.available_list,
            &mut self.processing_list,
        );

        // reset inverse_state
        self.batch_adder.reset();

        // previous points were processed except the last
        self.cur_points[0] = self.cur_points[self.cur_points_cnt as usize - 1];
        self.cur_points_cnt = 1;

        self.cur_batch_cnt = 0;
        self.collision_cnt = 0;

        // Process collision points
        // iterate over the unprocessed list
        //   - add to processing if no collision
        //   - add back to unprocessed if collision again
        let mut cur_node_idx = list::get_head(self.unprocessed_list);
        list::clear(&mut self.unprocessed_list);
        while cur_node_idx != 0x3FFF {
            // <bucket_id, next entry index> pair
            // cur_node_idx corresponding to a bucket_id in
            // collision_list_nodes and a point in cur_points+max_batch_cnt
            let cur_node = self.collision_list_nodes[cur_node_idx as usize];
            let bucket_id = list::get_entry_payload(cur_node);

            if !self.bitmap.test_and_set(bucket_id) {
                // print!("no collision\n");
                // if no collision get point from the collision area in cur_points
                let point = &self.cur_points[(self.max_batch_cnt + cur_node_idx) as usize];

                let acc = &self.buckets[bucket_id as usize];
                // space after max_batch_cnt for collisioned points
                let point_idx = self.max_batch_cnt + cur_node_idx;
                // print!("_process_batch: bucket_id {}, point index: {}, batch index {}\n", bucket_id, point_idx, self.cur_batch_cnt);
                self.batch_adder.batch_add_phase_one(acc, point, self.cur_batch_cnt as usize);

                // print!("bucket_id {}, point idx {}\n", bucket_id, cur_node_idx);

                self.cur_bkt_pnt_pairs[self.cur_batch_cnt as usize] =
                    ((bucket_id as u64) << 14) + point_idx as u64;
                self.cur_batch_cnt += 1;

                list::enqueue(
                    &mut self.collision_list_nodes,
                    &mut self.processing_list,
                    cur_node_idx,
                );
            } else {
                // print!("collision\n");
                self.collision_cnt += 1;
                list::enqueue(
                    &mut self.collision_list_nodes,
                    &mut self.unprocessed_list,
                    cur_node_idx,
                );
            }
            cur_node_idx = list::get_next(cur_node);
        }
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

    #[test]
    fn test_process_point_and_slices_deal_two_points() {
        let window_bits = 15u32;
        let mut bucket_msm = BucketMSM::<G1Affine>::new(30u32, window_bits, 128u32, 4096u32);
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
        let mut bucket_msm = BucketMSM::<G1Affine>::new(45u32, window_bits, 128u32, 4096u32);
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
}
