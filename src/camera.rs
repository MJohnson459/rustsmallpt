#![cfg_attr(feature = "unstable", feature(test))]

use vector_3d::Vec3d;
use ray::Ray;

pub struct Camera {
    pub ray: Ray,
    pub cx: Vec3d,
    pub cy: Vec3d
}

impl Camera {
    pub fn new(width: usize, height: usize, fov: f64) -> Camera {

        let ray: Ray = Ray{origin: Vec3d{x:50.0,y:50.0,z:295.6}, direction: Vec3d{x:0.0,y:-0.042612, z:-1.0}.normalise()};
        let cx: Vec3d = Vec3d{x:(width as f64)*fov/(height as f64),y:0.0,z:0.0}; // x direction increment
        let cy: Vec3d = (cx % ray.direction).normalise()*fov;                    // y direction increment

        Camera {
            ray: ray,
            cx: cx,
            cy: cy,
        }
    }
}
