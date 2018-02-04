pub struct Config {
    pub explicit_light_sampling: bool,
    pub height: usize,
    pub save_per_sample: bool,
    pub tent_filter: bool,
    pub width: usize,
    pub samples: usize,
    pub roulette_depth: i32,
}

impl Config {
    pub fn new(width: usize, height: usize, samples: usize) -> Config {
        Config {
            explicit_light_sampling: false,
            height: height,
            save_per_sample: true,
            tent_filter: true,
            width: width,
            samples: samples,
            roulette_depth: 5,
        }
    }
}
