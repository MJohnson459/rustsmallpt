#![cfg_attr(feature = "unstable", feature(test))]

extern crate rand;

use na::{Vector3};

use sphere::Sphere;
use material::ReflectType;
use ray::Ray;

#[derive(Debug, Deserialize)]
pub enum AvailableScenes {New1, New2}

pub struct Scene {
    pub name: String,
    pub spheres: Vec<Sphere>
}

impl Scene {
    pub fn new() -> Scene {
        let mut spheres = Vec::new();
        let wall_radius = 1e5;
        spheres.push(Sphere {radius: wall_radius, position: Vector3::<f64>::new(wall_radius + 1.0, 40.8, 81.6),    emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(0.75,0.25,0.25), reflection: ReflectType::DIFF}); // left
        spheres.push(Sphere {radius: wall_radius, position: Vector3::<f64>::new(-wall_radius + 99.0, 40.8, 81.6),  emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(0.25,0.25,0.75), reflection: ReflectType::DIFF}); // right
        spheres.push(Sphere {radius: wall_radius, position: Vector3::<f64>::new(50.0, 40.8, wall_radius),          emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(0.75,0.75,0.75), reflection: ReflectType::DIFF}); // back
        spheres.push(Sphere {radius: wall_radius, position: Vector3::<f64>::new(50.0, 40.8, -wall_radius + 170.0), emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(0.0,0.0,0.0), reflection: ReflectType::DIFF}); // front
        spheres.push(Sphere {radius: wall_radius, position: Vector3::<f64>::new(50.0, wall_radius, 81.6),          emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(0.75,0.75,0.75), reflection: ReflectType::DIFF}); // bottom
        spheres.push(Sphere {radius: wall_radius, position: Vector3::<f64>::new(50.0, -wall_radius + 81.6, 81.6),  emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(0.75,0.75,0.75), reflection: ReflectType::DIFF}); // Top
        spheres.push(Sphere {radius:16.5,         position: Vector3::<f64>::new(27.0, 16.5, 47.0),                 emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(1.0,1.0,1.0)*0.999, reflection: ReflectType::SPEC}); // mirror
        spheres.push(Sphere {radius:16.5,         position: Vector3::<f64>::new(73.0, 16.5, 78.0),                 emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(1.0,1.0,1.0)*0.999, reflection: ReflectType::REFR}); // glass
        spheres.push(Sphere {radius:600.0,        position: Vector3::<f64>::new(50.0, 681.6 - 0.27, 81.6),         emission: Vector3::<f64>::new(12.0,12.0,12.0), color: Vector3::<f64>::zeros(), reflection: ReflectType::DIFF}); // light

        Scene {
            name: String::from("New1"),
            spheres: spheres,
        }
    }

    pub fn new2() -> Scene {
        let mut spheres = Vec::new();
        let wall_radius = 1e5;
        let room_width = 100.0; // x
        let room_height = 80.0; // y
        let room_depth = 160.0; // z
        let camera_space = 100.0;
        let light_radius = 600.0;
        spheres.push(Sphere {radius: wall_radius,  position: Vector3::<f64>::new(wall_radius+1.0, room_height/2.0, room_depth/2.0),                      emission: Vector3::<f64>::zeros(),    color: Vector3::<f64>::new(0.8,0.1,0.1),    reflection: ReflectType::DIFF}); // left
        spheres.push(Sphere {radius: wall_radius,  position: Vector3::<f64>::new(-wall_radius+room_width-1.0, room_height/2.0, room_depth/2.0),          emission: Vector3::<f64>::zeros(),    color: Vector3::<f64>::new(0.1,0.3,0.70),   reflection: ReflectType::DIFF}); // right
        spheres.push(Sphere {radius: wall_radius,  position: Vector3::<f64>::new(room_width/2.0, room_height/2.0, wall_radius),                          emission: Vector3::<f64>::zeros(),    color: Vector3::<f64>::new(0.75,0.75,0.75), reflection: ReflectType::DIFF}); // back
        spheres.push(Sphere {radius: wall_radius,  position: Vector3::<f64>::new(room_width/2.0, room_height/2.0, -wall_radius+room_depth+camera_space), emission: Vector3::<f64>::zeros(),    color: Vector3::<f64>::new(0.0,0.0,0.0),    reflection: ReflectType::DIFF}); // front
        spheres.push(Sphere {radius: wall_radius,  position: Vector3::<f64>::new(room_width/2.0, wall_radius, room_depth/2.0),                           emission: Vector3::<f64>::zeros(),    color: Vector3::<f64>::new(0.75,0.75,0.75), reflection: ReflectType::DIFF}); // bottom
        spheres.push(Sphere {radius: wall_radius,  position: Vector3::<f64>::new(room_width/2.0, -wall_radius+room_height, room_depth/2.0),              emission: Vector3::<f64>::zeros(),    color: Vector3::<f64>::new(0.75,0.75,0.75), reflection: ReflectType::DIFF}); // Top
        spheres.push(Sphere {radius: light_radius, position: Vector3::<f64>::new(room_width/2.0, light_radius+room_height-0.27, room_depth/2.0),         emission: Vector3::<f64>::new(12.0,12.0,12.0), color: Vector3::<f64>::zeros(),    reflection: ReflectType::DIFF}); // light

        spheres.push(Sphere {radius: 16.5, position: Vector3::<f64>::new(27.0, 16.5, 47.0), emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(1.0,1.0,1.0)*0.95, reflection: ReflectType::SPEC}); // mirror
        spheres.push(Sphere {radius: 16.5, position: Vector3::<f64>::new(73.0, 30.0, 78.0), emission: Vector3::<f64>::zeros(), color: Vector3::<f64>::new(1.0,1.0,1.0)*0.95, reflection: ReflectType::REFR}); // glass

        Scene {
            name: String::from("New2"),
            spheres: spheres,
        }
    }

    pub fn intersect(&self, ray: &Ray) -> (bool, f64, usize) {
        let inf = 1e20;
        let mut closest_distance: f64 = inf;
        let mut id: usize = 0;

        for i in (0..self.spheres.len()).rev() {
            let distance = self.spheres[i].intersect(ray);
            if distance != 0.0 && distance < closest_distance {
                closest_distance = distance;
                id = i;
            }
        }
        (closest_distance < inf, closest_distance, id)
    }
}
