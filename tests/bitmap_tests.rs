#[cfg(test)]
mod bitmap_tests {
    use all_asserts::assert_false;
    use ark_msm::bitmap::Bitmap;

    fn setup() -> Bitmap {
        Bitmap::new(67)
    }

    #[test]
    fn test_set_when_only_once_return_false() {
        let mut bitmap = setup();
        assert_false!(bitmap.test_and_set(0));
        assert_false!(bitmap.test_and_set(33));
    }

    #[test]
    fn test_set_when_more_than_once_return_true() {
        let mut bitmap = setup();
        assert_false!(bitmap.test_and_set(0));
        assert!(bitmap.test_and_set(0));
    }

    #[test]
    fn test_clear() {
        let mut bitmap = setup();
        assert_false!(bitmap.test_and_set(0));
        assert!(bitmap.test_and_set(0));

        bitmap.clear();
        assert_false!(bitmap.test_and_set(0));
    }
}
