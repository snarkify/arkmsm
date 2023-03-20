#[cfg(test)]
mod bitmap_tests {
    use ark_msm::bitmap::Bitmap;

    fn setup() -> Bitmap {
        return Bitmap::new(67);
    }

    #[test]
    fn test_set_when_only_once_return_false() {
        let mut bitmap = setup();
        assert_eq!(bitmap.test_and_set(0), false);
        assert_eq!(bitmap.test_and_set(33), false);
    }

    #[test]
    fn test_set_when_more_than_once_return_true() {
        let mut bitmap = setup();
        assert_eq!(bitmap.test_and_set(0), false);
        assert_eq!(bitmap.test_and_set(0), true);
    }

    #[test]
    fn test_clear() {
        let mut bitmap = setup();
        assert_eq!(bitmap.test_and_set(0), false);
        assert_eq!(bitmap.test_and_set(0), true);

        bitmap.clear();
        assert_eq!(bitmap.test_and_set(0), false);
    }
}
