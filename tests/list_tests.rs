#[cfg(test)]
mod list_tests {
    use ark_msm::list::*;

    #[test]
    fn test_make_list() {
        let head = 1;
        let tail = 2;
        let list = make_list(head, tail);
        let head_result = get_head(list);
        let tail_result = get_tail(list);

        assert_eq!(head, head_result);
        assert_eq!(tail, tail_result);
    }

    #[test]
    fn test_enqueue() {
        let mut entries = vec![0; 10];
        let mut list = make_list(NIL, NIL);
        assert_eq!(get_head(list), NIL);

        enqueue(&mut entries, &mut list, 0);
        assert_eq!(get_head(list), 0);
        assert_eq!(get_tail(list), 0);

        enqueue(&mut entries, &mut list, 1);
        assert_eq!(get_head(list), 0);
        assert_eq!(get_tail(list), 1);

        enqueue(&mut entries, &mut list, 2);
        assert_eq!(get_head(list), 0);
        assert_eq!(get_tail(list), 2);
    }

    #[test]
    fn test_pop() {
        let mut entries = vec![0; 10];
        let mut list = make_list(NIL, NIL);

        enqueue(&mut entries, &mut list, 0);
        enqueue(&mut entries, &mut list, 1);
        enqueue(&mut entries, &mut list, 2);

        assert_eq!(pop(&mut entries, &mut list), 0);
        assert_eq!(pop(&mut entries, &mut list), 1);
        assert_eq!(pop(&mut entries, &mut list), 2);
        assert_eq!(pop(&mut entries, &mut list), NIL);
    }

    #[test]
    fn test_append_all_when_dest_is_empty() {
        let mut entries = vec![0; 10];
        let mut dest_list = make_list(NIL, NIL);
        let mut src_list = make_list(NIL, NIL);

        enqueue(&mut entries, &mut src_list, 0);
        enqueue(&mut entries, &mut src_list, 1);
        enqueue(&mut entries, &mut src_list, 2);
        enqueue(&mut entries, &mut src_list, 3);

        append_all(&mut entries, &mut dest_list, &mut src_list);

        assert_eq!(pop(&mut entries, &mut dest_list), 0);
        assert_eq!(pop(&mut entries, &mut dest_list), 1);
        assert_eq!(pop(&mut entries, &mut dest_list), 2);
        assert_eq!(pop(&mut entries, &mut dest_list), 3);
    }

    #[test]
    fn test_append_all_when_src_is_empty() {
        let mut entries = vec![0; 10];
        let mut src_list = make_list(NIL, NIL);
        let mut dest_list = make_list(NIL, NIL);

        enqueue(&mut entries, &mut dest_list, 0);
        enqueue(&mut entries, &mut dest_list, 1);
        enqueue(&mut entries, &mut dest_list, 2);
        enqueue(&mut entries, &mut dest_list, 3);

        append_all(&mut entries, &mut dest_list, &mut src_list);

        assert_eq!(pop(&mut entries, &mut dest_list), 0);
        assert_eq!(pop(&mut entries, &mut dest_list), 1);
        assert_eq!(pop(&mut entries, &mut dest_list), 2);
        assert_eq!(pop(&mut entries, &mut dest_list), 3);
    }

    #[test]
    fn test_append_all_when_both_are_nonempty() {
        let mut entries = vec![0; 10];
        let mut src_list = make_list(NIL, NIL);
        let mut dest_list = make_list(NIL, NIL);

        enqueue(&mut entries, &mut src_list, 0);
        enqueue(&mut entries, &mut src_list, 1);
        enqueue(&mut entries, &mut dest_list, 2);
        enqueue(&mut entries, &mut dest_list, 3);

        append_all(&mut entries, &mut dest_list, &mut src_list);

        assert_eq!(pop(&mut entries, &mut dest_list), 0);
        assert_eq!(pop(&mut entries, &mut dest_list), 1);
        assert_eq!(pop(&mut entries, &mut dest_list), 2);
        assert_eq!(pop(&mut entries, &mut dest_list), 3);
    }
}
