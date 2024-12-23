struct CIEDE2000 {}

impl CIEDE2000 {
    fn new() -> Self {
        CIEDE2000 {}
    }

    fn rgb565_to_rgb888(
        color: Rgb565,
    ) -> Rgb888 {
        Rgb888::new(
            color.r() << 3,
            color.g() << 2,
            color.b() << 3,
        )
    }

    fn rgb_to_xyz(
        color: Rgb888,
    ) -> (f32, f32, f32) {
        let r = color.r() as f32 / 255.0;
        let g = color.g() as f32 / 255.0;
        let b = color.b() as f32 / 255.0;

        let r = if r > 0.04045 {
            ((r + 0.055) / 1.055).powf(2.4)
        } else {
            r / 12.92
        };
        let g = if g > 0.04045 {
            ((g + 0.055) / 1.055).powf(2.4)
        } else {
            g / 12.92
        };
        let b = if b > 0.04045 {
            ((b + 0.055) / 1.055).powf(2.4)
        } else {
            b / 12.92
        };

        let x = r * 0.4124 + g * 0.3576 + b * 0.1805;
        let y = r * 0.2126 + g * 0.7152 + b * 0.0722;
        let z = r * 0.0193 + g * 0.1192 + b * 0.9505;

        (x, y, z)
    }

    fn xyz_to_lab(
        x: f32,
        y: f32,
        z: f32,
    ) -> (f32, f32, f32) {
        let x = x / 0.95047;
        let y = y / 1.00000;
        let z = z / 1.08883;

        let x = if x > 0.008856 {
            x.powf(1.0 / 3.0)
        } else {
            (7.787 * x) + (16.0 / 116.0)
        };
        let y = if y > 0.008856 {
            y.powf(1.0 / 3.0)
        } else {
            (7.787 * y) + (16.0 / 116.0)
        };
        let z = if z > 0.008856 {
            z.powf(1.0 / 3.0)
        } else {
            (7.787 * z) + (16.0 / 116.0)
        };

        let l = (116.0 * y) - 16.0;
        let a = 500.0 * (x - y);
        let b = 200.0 * (y - z);

        (l, a, b)
    }

    fn cie_de2000(
        lab1: (f32, f32, f32),
        lab2: (f32, f32, f32),
    ) -> f32 {
        let (l1, a1, b1) = lab1;
        let (l2, a2, b2) = lab2;

        let avg_lp = (l1 + l2) / 2.0;
        let c1 = (a1.powi(2) + b1.powi(2)).sqrt();
        let c2 = (a2.powi(2) + b2.powi(2)).sqrt();
        let avg_c = (c1 + c2) / 2.0;
        let g = 0.5 * (1.0 - (avg_c.powi(7)
            / (avg_c.powi(7) + 25.0f32.powi(7))).sqrt());

        let a1p = (1.0 + g) * a1;
        let a2p = (1.0 + g) * a2;
        let c1p = (a1p.powi(2) + b1.powi(2)).sqrt();
        let c2p = (a2p.powi(2) + b2.powi(2)).sqrt();

        let avg_cp = (c1p + c2p) / 2.0;
        let h1p = if b1 == 0.0 && a1p == 0.0 {
            0.0
        } else {
            b1.atan2(a1p).to_degrees() % 360.0
        };
        let h2p = if b2 == 0.0 && a2p == 0.0 {
            0.0
        } else {
            b2.atan2(a2p).to_degrees() % 360.0
        };

        let avg_hp = if (h1p - h2p).abs() > 180.0 {
            (h1p + h2p + 360.0) / 2.0
        } else {
            (h1p + h2p) / 2.0
        };

        let t = 1.0 - 0.17 * (avg_hp - 30.0).to_radians().cos()
            + 0.24 * (2.0 * avg_hp).to_radians().cos()
            + 0.32 * (3.0 * avg_hp + 6.0).to_radians().cos()
            - 0.20 * (4.0 * avg_hp - 63.0).to_radians().cos();

        let delta_hp = if h2p - h1p > 180.0 {
            h2p - h1p - 360.0
        } else if h2p - h1p < -180.0 {
            h2p - h1p + 360.0
        } else {
            h2p - h1p
        };

        let delta_lp = l2 - l1;
        let delta_cp = c2p - c1p;
        let delta_hp = 2.0 * (c1p * c2p).sqrt()
            * (delta_hp.to_radians() / 2.0).sin();

        let s_l = 1.0 + ((0.015 * (avg_lp - 50.0).powi(2))
            / ((20.0 + (avg_lp - 50.0).powi(2)).sqrt()));
        let s_c = 1.0 + 0.045 * avg_cp;
        let s_h = 1.0 + 0.015 * avg_cp * t;

        let delta_ro = 30.0 * (-((avg_hp - 275.0) / 25.0).powi(2)).exp();
        let r_c = 2.0 * (avg_cp.powi(7)
            / (avg_cp.powi(7) + 25.0f32.powi(7))).sqrt();
        let r_t = -r_c * (2.0 * delta_ro).to_radians().sin();

        ((delta_lp / s_l).powi(2)
            + (delta_cp / s_c).powi(2)
            + (delta_hp / s_h).powi(2)
            + r_t * (delta_cp / s_c) * (delta_hp / s_h))
            .sqrt()
    }

    fn ciede2000_distance(
        color1: Rgb565,
        color2: Rgb565,
    ) -> f32 {
        let rgb1 = Self::rgb565_to_rgb888(color1);
        let rgb2 = Self::rgb565_to_rgb888(color2);

        let (x1, y1, z1) = Self::rgb_to_xyz(rgb1);
        let (x2, y2, z2) = Self::rgb_to_xyz(rgb2);

        let lab1 = Self::xyz_to_lab(x1, y1, z1);
        let lab2 = Self::xyz_to_lab(x2, y2, z2);

        let distance = Self::cie_de2000(lab1, lab2);
        distance
    }
}
