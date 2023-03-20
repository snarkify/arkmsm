// TODO (ouyang):
// A new data structure should be created to deprecate lists in this file
// A structure named CollisionState should maintain two lists
//   - unprocessed: list of points observed in collision
//   - processing: list of points added to current batch that will be processed
//   and removed upon complete of processing

// Entries in a array are orgnized as multiple linked lists
//
// A list is reference by a <head, tail> pair, that is encoded in a u32
//    u32: (resv | 14 bit tail | 14 bit head)
//
// An entry in the list consists a payload and a next index
//     u32: (up to 28 bit payload | 14 bit next)
pub const NIL: u32 = 0x3FFF;
pub const MASK: u32 = NIL;

pub const fn make_list(head: u32, tail: u32) -> u32 {
    (tail << 14) | head
}

pub const fn get_head(list_ht: u32) -> u32 {
    list_ht & MASK
}

pub const fn get_tail(list_ht: u32) -> u32 {
    list_ht >> 14
}

pub const fn get_entry_payload(entry: u64) -> u32 {
    (entry >> 14) as u32
}

pub const fn get_next(entry: u64) -> u32 {
    (entry & MASK as u64) as u32
}

fn update_entry_next(entry: &mut u64, next: u32) {
    *entry = (*entry & 0xFFFFFFFFFFFFC000) | ((next & MASK) as u64);
}

pub fn is_empty(list_ht: u32) -> bool {
    list_ht == make_list(NIL, NIL)
}

pub fn clear(list: &mut u32) {
    *list = make_list(NIL, NIL)
}

// append entry[index] to list
pub fn enqueue(entries: &mut Vec<u64>, list: &mut u32, index: u32) {
    if is_empty(*list) {
        *list = make_list(index, index);
    } else {
        // update tail pointers
        let head = get_head(*list);
        let tail = get_tail(*list);
        *list = make_list(head, index);

        // append index to last entry
        update_entry_next(&mut entries[tail as usize], index);
    }
    entries[index as usize] |= NIL as u64;
}

// return index of a entry (with next pointer set)
pub fn pop(entries: &Vec<u64>, list: &mut u32) -> u32 {
    let list_old = *list;

    if !is_empty(*list) {
        // pop last entry in list
        if get_head(*list) == get_tail(*list) {
            clear(list);
        } else {
            let entry = entries[get_head(list_old) as usize];
            *list = make_list(get_next(entry), get_tail(list_old));
        }
    }

    get_head(list_old)
}

// insert src list to the head of the dest list
pub fn append_all(entries: &mut Vec<u64>, dest_list: &mut u32, src_list: &mut u32) {
    if is_empty(*src_list) {
        return;
    }

    if is_empty(*dest_list) {
        *dest_list = *src_list;
    } else {
        // point src tail to dest head
        update_entry_next(
            &mut entries[get_tail(*src_list) as usize],
            get_head(*dest_list),
        );

        *dest_list = make_list(get_head(*src_list), get_tail(*dest_list))
    }

    clear(src_list);
}
