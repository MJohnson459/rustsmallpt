use vector_3d::Vec3d;
use ray::Ray;
use material::ReflectType;

// ------
#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    pub radius: f64,
    pub position: Vec3d,
    pub emission: Vec3d,
    pub color: Vec3d,
    pub reflection: ReflectType
}

impl Sphere {

    // return distance 0.0 if nohit
    pub fn intersect(&self, ray: &Ray) -> f64 {
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        let op: Vec3d = self.position - ray.origin; // p is the sphere center (C)
        let eps: f64 = 1e-4;                        // eps is a small fudge factor
        let b: f64 = op.dot(&ray.direction);        // 1/2 b from quadratic eq. setup
        let mut det: f64 = b*b - op.dot(&op) + self.radius*self.radius; //(b^2-4ac)/4: a=1 because ray normalized
        if det < 0.0 {  // ray missed sphere
            return 0.0
        } else {
            det = det.sqrt();
        }

        // return smaller positive t
        let t1 = b - det;
        let t2 = b + det;
        if t1 > eps {
            return t1;
        } else if t2 > eps {
            return t2;
        } else {
            return 0.0;
        }
    }
}
