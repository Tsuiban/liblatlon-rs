use std::f64::consts::PI;

const R: f64 = 3443.9308855;

pub const NM_2_METERS: f64 = 1852.0;
pub const NM_2_FEET: f64 = 6076.12;

#[derive(Debug, Clone)]

pub struct DegreesDecimalMinutes {
    pub degrees: i8,
    pub minutes: f64,
}

pub struct DegreesMinutesSeconds {
    pub degrees: i8,
    pub minutes: u8,
    pub seconds: u8,
}

#[derive(Debug, Clone)]
pub struct Position {
    pub latitude: f64,
    pub longitude: f64,
}

impl From<(f64, f64)> for Position {
    fn from(pair: (f64, f64)) -> Self {
        Position {
            latitude: pair.0,
            longitude: pair.1,
        }
    }
}

impl From<(f32, f32)> for Position {
    fn from(pair: (f32, f32)) -> Self {
        Position {
            latitude: f64::from(pair.0),
            longitude: f64::from(pair.1),
        }
    }
}

impl Position {
    pub fn new() -> Position {
        Position {
            latitude: 0.0,
            longitude: 0.0,
        }
    }
    
    pub fn ne(&self, other: &Position) -> bool {
        self.latitude != other.latitude || self.longitude != other.longitude
    }

    pub fn from_degrees_decimal_minutes(
        latitude_degrees: i8,
        latitude_minutes: f64,
        longitude_degrees: i8,
        longitude_minutes: f64,
    ) -> Self {
        let latitude = f64::from(latitude_degrees)
            + (if latitude_degrees < 0 {
                -latitude_minutes
            } else {
                latitude_minutes
            }) / 60.0;
        let longitude = f64::from(longitude_degrees)
            + (if longitude_degrees < 0 {
                -longitude_minutes
            } else {
                longitude_minutes
            }) / 60.0;
        Position {
            latitude,
            longitude,
        }
    }

    pub fn from_hms(
        latitude_degrees: i8,
        latitude_minutes: u8,
        latitude_seconds: u8,
        longitude_degrees: i8,
        longitude_minutes: u8,
        longitude_seconds: u8,
    ) -> Self {
        let latitude = f64::from(latitude_degrees)
            + if latitude_degrees < 0 {
                -f64::from(latitude_minutes) / 60.0 - f64::from(latitude_seconds) / 3600.0
            } else {
                f64::from(latitude_minutes) / 60.0 + f64::from(latitude_seconds) / 3600.0
            };
        let longitude = f64::from(longitude_degrees)
            + if longitude_degrees < 0 {
                -f64::from(longitude_minutes) / 60.0 - f64::from(longitude_seconds) / 3600.0
            } else {
                f64::from(longitude_minutes) / 60.0 + f64::from(longitude_seconds) / 3600.0
            };
        Position {
            latitude,
            longitude,
        }
    }

    pub fn as_degrees_decimal_minutes(&self) -> (DegreesDecimalMinutes, DegreesDecimalMinutes) {
        let latitude = DegreesDecimalMinutes {
            degrees: self.latitude as i8,
            minutes: self.latitude.fract().abs() * 60.0,
        };
        let longitude = DegreesDecimalMinutes {
            degrees: self.longitude as i8,
            minutes: self.longitude.fract().abs() * 60.0,
        };
        (latitude, longitude)
    }

    pub fn as_degrees_minutes_seconds(&self) -> (DegreesMinutesSeconds, DegreesMinutesSeconds) {
        let minutes = self.latitude.fract().abs() * 60.0;
        let seconds = minutes.fract().abs() * 60.0;
        let latitude = DegreesMinutesSeconds {
            degrees: self.latitude as i8,
            minutes: minutes as u8,
            seconds: seconds as u8,
        };

        let minutes = self.longitude.fract().abs() * 60.0;
        let seconds = minutes.fract().abs() * 60.0;
        let longitude = DegreesMinutesSeconds {
            degrees: self.longitude as i8,
            minutes: minutes as u8,
            seconds: seconds as u8,
        };

        (latitude, longitude)
    }

    /// Calculates the distance to another position and returns the distance in nautical miles
    pub fn distance_to(&self, other: &Position) -> f64 {
        let phi_1 = self.latitude.to_radians();
        let phi_2 = other.latitude.to_radians();
        let delta_phi = phi_2 - phi_1;
        let delta_lambda = (other.longitude - self.longitude).to_radians();

        let delta_phi_sin = (delta_phi / 2.0).sin();
        let delta_lambda_sin = (delta_lambda / 2.0).sin();

        let a = (delta_phi_sin * delta_phi_sin)
            + phi_1.cos() * phi_2.cos() * (delta_lambda_sin * delta_lambda_sin);
        let c = 2.0 * (a.sqrt()).atan2((1.0 - a).sqrt());
        let d = R * c;
        d
    }

    /// Calculate the bearing to another position
    pub fn bearing_to(&self, other: &Position) -> f64 {
        let lambda_1 = self.longitude.to_radians();
        let lambda_2 = other.longitude.to_radians();
        let phi_1 = self.latitude.to_radians();
        let phi_2 = other.latitude.to_radians();

        let y = (lambda_2 - lambda_1).sin() * phi_2.cos();
        let x = phi_1.cos() * phi_2.sin() - phi_1.sin() * phi_2.cos() * (lambda_2 - lambda_1).cos();
        let theta = y.atan2(x);
        let brng = theta.to_degrees();
        if brng < 0.0 {
            brng + 360.0
        } else {
            brng
        }
    }

    /// find the midpoint between two positions.
    pub fn midpoint_between(&self, other: &Position) -> Position {
        let phi_1 = self.latitude.to_radians();
        let phi_2 = other.latitude.to_radians();
        let lambda_1 = self.longitude.to_radians();
        let lambda_2 = other.latitude.to_radians();

        let b_x = phi_2.cos() * (lambda_2 - lambda_1).cos();
        let b_y = phi_2.cos() * (lambda_2 - lambda_1).sin();
        let inner = phi_1.cos() + b_x;
        let phi_m = (phi_1.sin() + phi_2.sin()).atan2((inner * inner + b_y * b_y).sqrt());
        let lambda_m = lambda_1 + b_y.atan2(phi_1.cos() + b_x);

        Position {
            latitude: phi_m.to_degrees(),
            longitude: lambda_m.to_degrees(),
        }
    }

    /// Given a starting point, a course, and a distance, determine the ending position.

    pub fn destination(&self, bearing: f64, distance: f64) -> Position {
        let phi_1 = self.latitude.to_radians();
        let lambda_1 = self.longitude.to_radians();
        let theta = bearing.to_radians();
        let omega = R / distance;

        let phi_2 = phi_1.sin() * theta.cos() + phi_1.cos() * theta.sin() * omega.cos();
        let lambda_2 = lambda_1
            + (omega.sin() * theta.sin() * phi_1.cos())
                .atan2(theta.cos() - phi_1.sin() * phi_2.sin());
        Position {
            latitude: phi_2.to_degrees(),
            longitude: lambda_2.to_degrees(),
        }
    }

    /// Given two starting points and courses from those starting points, determine where they will
    /// intersect.
    pub fn intersection(&self, course_1: f64, other: &Position, course_2: f64) -> Position {
        let phi_1 = self.latitude.to_radians();
        let phi_2 = other.latitude.to_radians();
        let lambda_1 = self.longitude.to_radians();
        let lambda_2 = other.longitude.to_radians();
        let theta_13 = course_1.to_radians();
        let theta_23 = course_2.to_radians();

        let sin_delta_phi = ((phi_2 - phi_1) / 2.0).sin();
        let sin_delta_lambda = ((lambda_2 - lambda_1) / 2.0).sin();

        let sigma_12 = (sin_delta_phi * sin_delta_phi
            + phi_1.cos() * phi_2.cos() * sin_delta_lambda * sin_delta_lambda)
            .sqrt()
            .asin()
            * 2.0;
        let theta_a =
            ((phi_2.sin() - phi_1.sin() * sigma_12.cos()) / (sigma_12.sin() * phi_1.cos())).acos();
        let theta_b =
            ((phi_1.sin() - phi_2.sin() * sigma_12.cos()) / (sigma_12.sin() * phi_2.cos())).acos();

        let theta_12: f64;
        let theta_21: f64;

        if (lambda_2 - lambda_1).sin() > 0.0 {
            theta_12 = theta_a;
            theta_21 = 2.0 * PI - theta_b;
        } else {
            theta_12 = 2.0 * PI - theta_a;
            theta_21 = theta_b;
        }

        let alpha_1 = theta_13 - theta_12;
        let alpha_2 = theta_21 - theta_23;

        let alpha_3 = (-(alpha_1.cos()) * alpha_2.cos()
            + alpha_1.sin() * alpha_2.sin() * sigma_12.cos())
        .acos();

        let sigma_13 = (sigma_12.sin() * alpha_1.sin() * alpha_2.sin())
            .atan2(alpha_2.cos() + alpha_1.cos() * alpha_3.cos());
        let phi_3 =
            (phi_1.sin() * sigma_13.cos() + phi_1.cos() * sigma_13.sin() * theta_13.cos()).asin();
        let delta_lambda_13 = (theta_13.sin() * sigma_13.sin() * phi_1.cos())
            .atan2(sigma_13.cos() - phi_1.sin() * phi_3.sin());
        let lambda_3 = lambda_1 + delta_lambda_13;

        Position {
            latitude: phi_3.to_degrees(),
            longitude: lambda_3.to_degrees(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CLOSE_TOLERANCE: f64 = 0.05;

    #[test]
    fn test01() {
        let position = Position::new();
        assert_eq!(position.latitude, 0.0);
        assert_eq!(position.longitude, 0.0);
    }

    #[test]
    fn test02() {
        let position = Position::from_degrees_decimal_minutes(1, 30.0, -120, 45.0);
        println!("{:?}", position);
        assert_eq!(position.latitude, 1.5);
        assert_eq!(position.longitude, -120.75);
    }

    #[test]
    fn test03() {
        let position = Position::from_hms(123, 45, 15, -100, 30, 55);
        println!("{:?}", position);
        assert!((position.latitude - 123.75 - 15.0 / 3600.0).abs() < 0.1);
        assert!((position.longitude - (-100.5 + 55.0 / 3600.0)).abs() < 0.1);
    }

    #[test]
    fn test04() {
        let start_position = Position::from((0.0, 0.0));
        let end_position = Position::from((1.0, 0.0));
        assert!(start_position.distance_to(&end_position) - 60.0 < 0.11);

        let end_position = Position::from((0.0, 1.0));
        assert!(start_position.distance_to(&end_position) - 60.0 < 0.11);

        let start_position = Position::from((1.0, 0.0));
        let distance = start_position.distance_to(&end_position);
        let check = 2.0_f64.sqrt() * 60.0;
        assert!((distance - check).abs() < 0.01 * check);
    }

    #[test]
    fn test05() {
        let start_position = Position::from((0.0, 0.0));
        let end_position = Position::from((1.0, 0.0));
        let bearing = start_position.bearing_to(&end_position);

        assert!(bearing.abs() < CLOSE_TOLERANCE);

        let bearing = start_position.bearing_to(&Position::from((0.0, 1.0)));
        assert!((bearing - 90.0).abs() < CLOSE_TOLERANCE);

        let bearing = start_position.bearing_to(&Position::from((-1.0, 0.0)));
        assert!((bearing - 180.0).abs() < CLOSE_TOLERANCE);

        let bearing = start_position.bearing_to(&Position::from((0.0, -1.0)));
        assert!((bearing - 270.0).abs() < CLOSE_TOLERANCE);

        let bearing = start_position.bearing_to(&Position::from((1.0, 1.0)));
        assert!((bearing - 45.0).abs() < CLOSE_TOLERANCE);

        let bearing = start_position.bearing_to(&Position::from((-1.0, 1.0)));
        assert!((bearing - 135.0).abs() < CLOSE_TOLERANCE);

        let bearing = start_position.bearing_to(&Position::from((1.0, -1.0)));
        assert!((bearing - 315.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06a() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((49.0, 124.0)));
        println!("bearing is {bearing}");
        assert!((bearing - 270.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06b() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((0.0, 1.0)));
        assert!((bearing - 90.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06c() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((-1.0, 0.0)));
        assert!((bearing - 180.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06d() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((0.0, -1.0)));
        assert!((bearing - 270.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06e() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((1.0, 1.0)));
        assert!((bearing - 45.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06f() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((-1.0, 1.0)));
        assert!((bearing - 135.0).abs() < CLOSE_TOLERANCE);
    }

    #[test]
    fn test06g() {
        let start_position = Position::from((49.0, 123.0));

        let bearing = start_position.bearing_to(&Position::from((1.0, -1.0)));
        assert!((bearing - 315.0).abs() < CLOSE_TOLERANCE);
    }
}
