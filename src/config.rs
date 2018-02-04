use std::path::PathBuf;
use scene::{Scene, AvailableScenes};


pub struct Config {
    pub explicit_light_sampling: bool,
    pub height: usize,
    pub path: PathBuf,
    pub roulette_depth: i32,
    pub samples: usize,
    pub save_per_sample: bool,
    pub scene: Scene,
    pub tent_filter: bool,
    pub width: usize,
}

impl Config {
    pub fn new(width: usize, height: usize, samples: usize, scene_name: Option<AvailableScenes>) -> Config {

        let scene =
            match scene_name {
                Some(x) => {
                    Scene::from_available(x)
                },
                None => {
                    println!("Using default scene Floating");
                    Scene::from_available(AvailableScenes::Floating)
                }
            };

        let image_name = format!("{}_{}_{}_{}.png", scene.name, width, height, samples);
        let path = PathBuf::from(&image_name);

        Config {
            explicit_light_sampling: false,
            height: height,
            path: path,
            roulette_depth: 5,
            samples: samples,
            save_per_sample: true,
            scene: scene,
            tent_filter: true,
            width: width,
        }
    }
}
