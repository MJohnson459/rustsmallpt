#![cfg_attr(feature = "unstable", feature(test))]

// ----------
#[derive(Copy, Clone, Debug)]
pub enum ReflectType {
    DIFF,
    SPEC,
    REFR
}
