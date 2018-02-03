extern crate rand;

use vector_3d::Vec3d;
use sphere::Sphere;
use material::ReflectType;
use ray::Ray;

#[derive(Debug, Deserialize)]
pub enum AvailableScenes {New1, New2}

pub struct Scene {
    pub name: String,
    pub spheres: Vec<Sphere>,
}

impl Scene {
    pub fn new() -> Scene {
        let mut spheres = Vec::new();
        let wall_radius = 1e5;
        spheres.push(Sphere {radius: wall_radius, position: Vec3d{x: wall_radius + 1.0, y: 40.8, z: 81.6},    emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.25,z:0.25}, reflection: ReflectType::DIFF}); // left
        spheres.push(Sphere {radius: wall_radius, position: Vec3d{x: -wall_radius + 99.0, y: 40.8, z: 81.6},  emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.25,y:0.25,z:0.75}, reflection: ReflectType::DIFF}); // right
        spheres.push(Sphere {radius: wall_radius, position: Vec3d{x: 50.0, y: 40.8, z: wall_radius},          emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // back
        spheres.push(Sphere {radius: wall_radius, position: Vec3d{x: 50.0, y: 40.8, z: -wall_radius + 170.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // front
        spheres.push(Sphere {radius: wall_radius, position: Vec3d{x: 50.0, y: wall_radius, z: 81.6},          emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // bottom
        spheres.push(Sphere {radius: wall_radius, position: Vec3d{x: 50.0, y: -wall_radius + 81.6, z: 81.6},  emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // Top
        spheres.push(Sphere {radius: 16.5,        position: Vec3d{x: 27.0, y: 16.5, z: 47.0},                 emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.999, reflection: ReflectType::SPEC}); // mirror
        spheres.push(Sphere {radius: 16.5,        position: Vec3d{x: 73.0, y: 16.5, z: 78.0},                 emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.999, reflection: ReflectType::REFR}); // glass
        spheres.push(Sphere {radius: 600.0,       position: Vec3d{x: 50.0, y: 681.6 - 0.27, z: 81.6},         emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light

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
        spheres.push(Sphere {radius: wall_radius,  position: Vec3d{x: wall_radius+1.0, y: room_height/2.0, z: room_depth/2.0},                      emission: Vec3d{x:0.0,y:0.0,z:0.0},    color: Vec3d{x:0.8,y:0.1,z:0.1},    reflection: ReflectType::DIFF}); // left
        spheres.push(Sphere {radius: wall_radius,  position: Vec3d{x: -wall_radius+room_width-1.0, y: room_height/2.0, z: room_depth/2.0},          emission: Vec3d{x:0.0,y:0.0,z:0.0},    color: Vec3d{x:0.1,y:0.3,z:0.70},   reflection: ReflectType::DIFF}); // right
        spheres.push(Sphere {radius: wall_radius,  position: Vec3d{x: room_width/2.0, y: room_height/2.0, z: wall_radius},                          emission: Vec3d{x:0.0,y:0.0,z:0.0},    color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // back
        spheres.push(Sphere {radius: wall_radius,  position: Vec3d{x: room_width/2.0, y: room_height/2.0, z: -wall_radius+room_depth+camera_space}, emission: Vec3d{x:0.0,y:0.0,z:0.0},    color: Vec3d{x:0.0,y:0.0,z:0.0},    reflection: ReflectType::DIFF}); // front
        spheres.push(Sphere {radius: wall_radius,  position: Vec3d{x: room_width/2.0, y: wall_radius, z: room_depth/2.0},                           emission: Vec3d{x:0.0,y:0.0,z:0.0},    color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // bottom
        spheres.push(Sphere {radius: wall_radius,  position: Vec3d{x: room_width/2.0, y: -wall_radius+room_height, z: room_depth/2.0},              emission: Vec3d{x:0.0,y:0.0,z:0.0},    color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // Top
        spheres.push(Sphere {radius: light_radius, position: Vec3d{x: room_width/2.0, y: light_radius+room_height-0.27, z: room_depth/2.0},         emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0},    reflection: ReflectType::DIFF}); // light

        spheres.push(Sphere {radius: 16.5, position: Vec3d{x: 27.0, y: 16.5, z: 47.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.95, reflection: ReflectType::SPEC}); // mirror
        spheres.push(Sphere {radius: 16.5, position: Vec3d{x: 73.0, y: 30.0, z: 78.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.95, reflection: ReflectType::REFR}); // glass

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
