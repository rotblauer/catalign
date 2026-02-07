//! Track management and rendering

use crate::cali::TileMetrics;

/// Manages multiple tracks and their visibility
pub struct TrackManager {
    pub show_identity: bool,
    pub show_coverage: bool,
    pub show_quality: bool,
    pub show_gap_rate: bool,
    pub show_energy: bool,
}

impl TrackManager {
    pub fn new() -> Self {
        Self {
            show_identity: true,
            show_coverage: true,
            show_quality: true,
            show_gap_rate: true,
            show_energy: false,
        }
    }
}

impl Default for TrackManager {
    fn default() -> Self {
        Self::new()
    }
}
