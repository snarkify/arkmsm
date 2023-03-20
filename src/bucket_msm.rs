use ark_bls12_381::G1Affine;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::Field;
use ark_std::{One, Zero};

use crate::{bitmap::Bitmap, list, types::G1Projective};

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

    // inverse state
    inverse_state: <G1Affine as AffineCurve>::BaseField,
    inverses: Vec<<G1Affine as AffineCurve>::BaseField>, // size max_batch_cnt

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

            inverse_state: <G1Affine as AffineCurve>::BaseField::one(),
            inverses: vec![<G1Affine as AffineCurve>::BaseField::one(); batch_size as usize],

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
            BucketMSM::<G1Affine>::batch_add_phase_one(
                acc,
                point,
                self.cur_batch_cnt as usize,
                &mut self.inverse_state,
                &mut self.inverses,
            );

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
        self.inverse_state = self.inverse_state.inverse().unwrap();
        for i in (0..self.cur_batch_cnt).rev() {
            let bucket_point = self.cur_bkt_pnt_pairs[i as usize];
            let point = &self.cur_points[(bucket_point & 0x3FFF) as usize];
            let acc = &mut self.buckets[(bucket_point >> 14) as usize];
            // print!("phase_two: bucket_id: {}, pair index: {}, point_idx {}\n", bucket_point >> 14, i, bucket_point & 0x3FFF);
            // print!{"phase_two acc: {}\n", acc};
            // print!{"phase_two point: {}\n", point};
            BucketMSM::<G1Affine>::batch_add_phase_two(
                acc,
                point,
                i as usize,
                &mut self.inverse_state,
                &mut self.inverses,
            );
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
        self.inverse_state.set_one();

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
                BucketMSM::<G1Affine>::batch_add_phase_one(
                    acc,
                    point,
                    self.cur_batch_cnt as usize,
                    &mut self.inverse_state,
                    &mut self.inverses,
                );

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

    // Two-pass batch affine addition
    //   - 1st pass calculates from left to right
    //      - state: accumulated product of deltaX
    //      - inverses[]: accumulated product left to a point
    //   - inverse state
    //   - 2nd pass calculates from right to left
    //      - slope s and ss from state
    //      - state =  state * deltaX
    //      - addition result acc
    fn batch_add_phase_one(
        p: &G1Affine,
        q: &G1Affine,
        idx: usize,
        inverse_state: &mut <G1Affine as AffineCurve>::BaseField,
        inverses: &mut Vec<<G1Affine as AffineCurve>::BaseField>,
    ) {
        if p.is_zero() | q.is_zero() {
            return;
        }

        let mut delta_x = q.x - p.x;
        if delta_x.is_zero() {
            let delta_y = q.y - p.y;
            if !delta_y.is_zero() {
                // p = -q, return
                return;
            }

            // if P == Q
            // if delta_x is zero, we need to invert 2y
            delta_x = q.y + q.y;
        }

        if inverse_state.is_zero() {
            inverses[idx].set_one();
            *inverse_state = delta_x;
        } else {
            inverses[idx] = *inverse_state;
            *inverse_state *= delta_x
        }
    }

    // should call compute inverse of state->inverse_state between phase_one and phase_two
    fn batch_add_phase_two(
        p: &mut G1Affine,
        q: &G1Affine,
        idx: usize,
        inverse_state: &mut <G1Affine as AffineCurve>::BaseField,
        inverses: &Vec<<G1Affine as AffineCurve>::BaseField>,
    ) {
        if p.is_zero() | q.is_zero() {
            if !q.is_zero() {
                *p = q.clone();
            }
            return;
        }

        let mut _inverse = inverses[idx];
        _inverse *= *inverse_state;

        let mut delta_x = q.x - p.x;
        let mut delta_y = q.y - p.y;

        if delta_x.is_zero() {
            if !delta_y.is_zero() {
                // p = -q, result should be pt at infinity
                p.set_zero();
                return;
            }
            // Otherwise, p = q, and it's point doubling
            // Processing is almost the same, except s=3*affine.x^2 / 2*affine.y

            // set delta_y = 3*q.x^2
            delta_y = q.x.square();
            delta_y = delta_y + delta_y + delta_y;

            delta_x = q.y.double();
        }

        // get the state ready for the next iteration
        *inverse_state *= delta_x;

        let s = delta_y * _inverse;
        let ss = s * s;
        p.x = ss - q.x - p.x;
        delta_x = q.x - p.x;
        p.y = s * delta_x;
        p.y = p.y - q.y;
    }

    pub fn msm_reduce(&mut self) -> G1Projective {
        let window_starts: Vec<_> = (0..self.num_windows).collect();

        let window_sums: Vec<G1Projective> = ark_std::cfg_into_iter!(window_starts)
            .map(|w_start| {
                let bucket_start = (w_start << self.bucket_bits) as usize;
                let bucket_end = (bucket_start + (1 << self.bucket_bits)) as usize;
                self.inner_window_reduce(bucket_start, bucket_end)
            })
            .collect();

        return self.intra_window_reduce(&window_sums);
    }

    fn inner_window_reduce(&mut self, start: usize, end: usize) -> G1Projective {
        let mut running_sum = G1Projective::zero();
        let mut res = G1Projective::zero();
        self.buckets[start..end].iter().rev().for_each(|b| {
            running_sum.add_assign_mixed(b);
            res += &running_sum;
        });
        return res;
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

#[cfg(test)]
mod batch_add_tests {
    use super::*;
    use ark_ec::ProjectiveCurve;
    use ark_std::UniformRand;
    use std::ops::Add;

    #[test]
    fn test_phase_one_zero_or_neg() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];
        BucketMSM::<G1Affine>::batch_add_phase_one(
            &G1Affine::zero(),
            &G1Affine::zero(),
            0,
            &mut inverse_state,
            &mut inverses,
        );

        let mut rng = ark_std::test_rng();
        let p = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p_affine = G1Affine::from(p);
        let mut neg_p_affine = p_affine.clone();
        neg_p_affine.y = -neg_p_affine.y;

        BucketMSM::<G1Affine>::batch_add_phase_one(
            &p_affine,
            &neg_p_affine,
            0,
            &mut inverse_state,
            &mut inverses,
        );
    }

    #[test]
    fn test_phase_one_p_add_p() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];
        let mut rng = ark_std::test_rng();
        let prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p = G1Affine::from(prj);
        let acc = p.clone();

        BucketMSM::<G1Affine>::batch_add_phase_one(&acc, &p, 0, &mut inverse_state, &mut inverses);
        assert_eq!(inverses[0].is_one(), true);
        assert_eq!(inverse_state, p.y + p.y);
    }

    #[test]
    fn test_phase_one_p_add_q() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];
        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let q_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p = G1Affine::from(p_prj);
        let q = G1Affine::from(q_prj);

        BucketMSM::<G1Affine>::batch_add_phase_one(&p, &q, 0, &mut inverse_state, &mut inverses);
        assert_eq!(inverses[0].is_one(), true);
        assert_eq!(inverse_state, q.x - p.x);
    }

    #[test]
    fn test_phase_one_p_add_q_twice() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];
        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let q_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p = G1Affine::from(p_prj);
        let q = G1Affine::from(q_prj);

        BucketMSM::<G1Affine>::batch_add_phase_one(&p, &q, 0, &mut inverse_state, &mut inverses);
        BucketMSM::<G1Affine>::batch_add_phase_one(&p, &q, 0, &mut inverse_state, &mut inverses);
        assert_eq!(inverses[0], q.x - p.x);
        assert_eq!(inverse_state, (q.x - p.x) * (q.x - p.x));
    }

    #[test]
    fn test_phase_two_zero_add_p() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];
        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let p = G1Affine::from(p_prj);
        let mut acc = G1Affine::zero();
        BucketMSM::<G1Affine>::batch_add_phase_two(
            &mut acc,
            &p,
            0,
            &mut inverse_state,
            &mut inverses,
        );
        assert_eq!(acc, p);
    }

    #[test]
    fn test_phase_two_p_add_neg() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];

        let mut rng = ark_std::test_rng();
        let p_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let mut acc = G1Affine::from(p_prj);
        let mut p = acc.clone();
        p.y = -p.y;

        BucketMSM::<G1Affine>::batch_add_phase_two(
            &mut acc,
            &p,
            0,
            &mut inverse_state,
            &mut inverses,
        );
        assert_eq!(acc, G1Affine::zero());
    }

    #[test]
    fn test_phase_two_p_add_q() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];
        let mut bucket_msm = BucketMSM::<G1Affine>::new(255u32, 15u32, 128u32, 4096u32);

        let mut rng = ark_std::test_rng();
        let acc_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let mut acc = G1Affine::from(acc_prj);
        let mut p = acc.clone();
        p.x = p.x + p.x;

        bucket_msm.inverses[0] = (p.x - acc.x).inverse().unwrap();
        BucketMSM::<G1Affine>::batch_add_phase_two(
            &mut acc,
            &p,
            0,
            &mut inverse_state,
            &mut inverses,
        );
        assert_eq!(acc, G1Affine::from(acc_prj.add_mixed(&p)));
    }

    #[test]
    fn test_phase_two_p_add_p() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 4 as usize];

        let mut rng = ark_std::test_rng();
        let acc_prj = <G1Affine as AffineCurve>::Projective::rand(&mut rng);
        let mut acc = G1Affine::from(acc_prj);
        let p = acc.clone();

        inverses[0] = (p.y + p.y).inverse().unwrap();
        BucketMSM::<G1Affine>::batch_add_phase_two(
            &mut acc,
            &p,
            0,
            &mut inverse_state,
            &mut inverses,
        );
        assert_eq!(acc, G1Affine::from(acc_prj).add(p));
    }

    #[test]
    fn test_batch_add() {
        let mut inverse_state = <G1Affine as AffineCurve>::BaseField::one();
        let mut inverses = vec![<G1Affine as AffineCurve>::BaseField::one(); 10 as usize];

        let mut rng = ark_std::test_rng();
        let mut buckets: Vec<G1Affine> = (0..10)
            .map(|_| G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(&mut rng)))
            .collect();
        let points: Vec<G1Affine> = (0..10)
            .map(|_| G1Affine::from(<G1Affine as AffineCurve>::Projective::rand(&mut rng)))
            .collect();

        let tmp = buckets.clone();
        for i in 0..10 {
            BucketMSM::<G1Affine>::batch_add_phase_one(
                &buckets[i],
                &points[i],
                i,
                &mut inverse_state,
                &mut inverses,
            );
        }
        inverse_state = inverse_state.inverse().unwrap();
        for i in (0..10).rev() {
            BucketMSM::<G1Affine>::batch_add_phase_two(
                &mut buckets[i],
                &points[i],
                i,
                &mut inverse_state,
                &mut inverses,
            );
        }

        for i in 0..10 {
            assert_eq!(buckets[i], tmp[i].add(points[i]));
        }
    }
}
