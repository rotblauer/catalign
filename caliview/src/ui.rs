//! UI components for track rendering

use egui::{Color32, Pos2, Rect, Stroke, Ui};
use crate::cali::TileMetrics;

/// Draw the identity track
pub fn draw_identity_track(
    ui: &mut Ui,
    tiles: &[TileMetrics],
    view_range: (u64, u64),
    height: f32,
) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), height),
        egui::Sense::hover(),
    );

    let rect = response.rect;
    let range_width = (view_range.1 - view_range.0) as f32;

    // Background
    painter.rect_filled(rect, 0.0, Color32::from_gray(20));

    // Draw tiles
    for tile in tiles {
        if tile.coverage == 0.0 {
            continue;
        }

        let x_start = rect.left() + (tile.start as f32 - view_range.0 as f32) / range_width * rect.width();
        let x_end = rect.left() + (tile.end as f32 - view_range.0 as f32) / range_width * rect.width();

        // Color based on identity (green = good, red = bad)
        let color = identity_to_color(tile.identity);

        let bar_height = tile.identity * (height - 20.0);
        let bar_rect = Rect::from_min_max(
            Pos2::new(x_start, rect.bottom() - 10.0 - bar_height),
            Pos2::new(x_end.max(x_start + 1.0), rect.bottom() - 10.0),
        );

        painter.rect_filled(bar_rect, 0.0, color);
    }

    // Y-axis labels
    painter.text(
        Pos2::new(rect.left() + 5.0, rect.top() + 5.0),
        egui::Align2::LEFT_TOP,
        "100%",
        egui::FontId::proportional(10.0),
        Color32::WHITE,
    );
    painter.text(
        Pos2::new(rect.left() + 5.0, rect.bottom() - 15.0),
        egui::Align2::LEFT_BOTTOM,
        "0%",
        egui::FontId::proportional(10.0),
        Color32::WHITE,
    );
}

/// Draw the coverage track
pub fn draw_coverage_track(
    ui: &mut Ui,
    tiles: &[TileMetrics],
    view_range: (u64, u64),
    height: f32,
) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), height),
        egui::Sense::hover(),
    );

    let rect = response.rect;
    let range_width = (view_range.1 - view_range.0) as f32;

    painter.rect_filled(rect, 0.0, Color32::from_gray(20));

    for tile in tiles {
        let x_start = rect.left() + (tile.start as f32 - view_range.0 as f32) / range_width * rect.width();
        let x_end = rect.left() + (tile.end as f32 - view_range.0 as f32) / range_width * rect.width();

        // Blue for coverage
        let intensity = (tile.coverage * 255.0) as u8;
        let color = Color32::from_rgb(50, 100, intensity);

        let bar_height = tile.coverage * (height - 20.0);
        let bar_rect = Rect::from_min_max(
            Pos2::new(x_start, rect.bottom() - 10.0 - bar_height),
            Pos2::new(x_end.max(x_start + 1.0), rect.bottom() - 10.0),
        );

        painter.rect_filled(bar_rect, 0.0, color);
    }

    // Labels
    painter.text(
        Pos2::new(rect.left() + 5.0, rect.top() + 5.0),
        egui::Align2::LEFT_TOP,
        "100%",
        egui::FontId::proportional(10.0),
        Color32::WHITE,
    );
}

/// Draw the quality score track
pub fn draw_quality_track(
    ui: &mut Ui,
    tiles: &[TileMetrics],
    view_range: (u64, u64),
    height: f32,
) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), height),
        egui::Sense::hover(),
    );

    let rect = response.rect;
    let range_width = (view_range.1 - view_range.0) as f32;

    painter.rect_filled(rect, 0.0, Color32::from_gray(20));

    for tile in tiles {
        if tile.coverage == 0.0 {
            continue;
        }

        let x_start = rect.left() + (tile.start as f32 - view_range.0 as f32) / range_width * rect.width();
        let x_end = rect.left() + (tile.end as f32 - view_range.0 as f32) / range_width * rect.width();

        // Purple gradient for quality
        let q = tile.quality_score / 100.0;
        let color = Color32::from_rgb(
            ((1.0 - q) * 200.0) as u8,
            (q * 100.0) as u8,
            (q * 200.0) as u8,
        );

        let bar_height = q * (height - 20.0);
        let bar_rect = Rect::from_min_max(
            Pos2::new(x_start, rect.bottom() - 10.0 - bar_height),
            Pos2::new(x_end.max(x_start + 1.0), rect.bottom() - 10.0),
        );

        painter.rect_filled(bar_rect, 0.0, color);
    }

    painter.text(
        Pos2::new(rect.left() + 5.0, rect.top() + 5.0),
        egui::Align2::LEFT_TOP,
        "Q100",
        egui::FontId::proportional(10.0),
        Color32::WHITE,
    );
}

/// Draw the gap rate track
pub fn draw_gap_track(
    ui: &mut Ui,
    tiles: &[TileMetrics],
    view_range: (u64, u64),
    height: f32,
) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), height),
        egui::Sense::hover(),
    );

    let rect = response.rect;
    let range_width = (view_range.1 - view_range.0) as f32;

    painter.rect_filled(rect, 0.0, Color32::from_gray(20));

    for tile in tiles {
        if tile.coverage == 0.0 {
            continue;
        }

        let x_start = rect.left() + (tile.start as f32 - view_range.0 as f32) / range_width * rect.width();
        let x_end = rect.left() + (tile.end as f32 - view_range.0 as f32) / range_width * rect.width();

        // Orange/red for gaps
        let gap = tile.gap_rate.min(1.0);
        let color = Color32::from_rgb(
            (200.0 + gap * 55.0) as u8,
            ((1.0 - gap) * 150.0) as u8,
            50,
        );

        let bar_height = gap * (height - 20.0);
        let bar_rect = Rect::from_min_max(
            Pos2::new(x_start, rect.bottom() - 10.0 - bar_height),
            Pos2::new(x_end.max(x_start + 1.0), rect.bottom() - 10.0),
        );

        painter.rect_filled(bar_rect, 0.0, color);
    }

    painter.text(
        Pos2::new(rect.left() + 5.0, rect.top() + 5.0),
        egui::Align2::LEFT_TOP,
        "Gap%",
        egui::FontId::proportional(10.0),
        Color32::WHITE,
    );
}

/// Convert identity to color (green = good, yellow = medium, red = bad)
fn identity_to_color(identity: f32) -> Color32 {
    if identity >= 0.99 {
        Color32::from_rgb(100, 200, 100)  // Green
    } else if identity >= 0.95 {
        Color32::from_rgb(150, 200, 80)   // Yellow-green
    } else if identity >= 0.90 {
        Color32::from_rgb(200, 200, 50)   // Yellow
    } else if identity >= 0.80 {
        Color32::from_rgb(220, 150, 50)   // Orange
    } else {
        Color32::from_rgb(200, 80, 80)    // Red
    }
}

/// Get color for a nucleotide base
fn base_to_color(base: char) -> Color32 {
    match base {
        'A' | 'a' => Color32::from_rgb(100, 180, 100),  // Green
        'C' | 'c' => Color32::from_rgb(100, 100, 200),  // Blue
        'G' | 'g' => Color32::from_rgb(200, 180, 100),  // Yellow/Gold
        'T' | 't' => Color32::from_rgb(200, 100, 100),  // Red
        'N' | 'n' => Color32::from_gray(100),           // Gray
        '-' => Color32::from_rgb(150, 50, 150),         // Purple for gaps
        _ => Color32::from_gray(80),
    }
}

/// Draw the sequence track with ACTG bases and indel highlighting
/// Only shown when zoomed in enough to see individual bases
pub fn draw_sequence_track(
    ui: &mut Ui,
    tiles: &[TileMetrics],
    view_range: (u64, u64),
    height: f32,
    reference_seq: Option<&str>,
    query_seq: Option<&str>,
) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), height),
        egui::Sense::hover(),
    );

    let rect = response.rect;
    let range_width = view_range.1 - view_range.0;
    let bases_per_pixel = range_width as f32 / rect.width();

    // Background
    painter.rect_filled(rect, 0.0, Color32::from_gray(25));

    // Only draw individual bases if zoomed in enough
    if bases_per_pixel <= 1.0 {
        // We can show individual bases
        let font_size = (rect.width() / range_width as f32).min(14.0).max(8.0);
        let row_height = height / 3.0;

        // Draw reference sequence (top row)
        if let Some(ref_seq) = reference_seq {
            let start_idx = view_range.0 as usize;
            let end_idx = (view_range.1 as usize).min(ref_seq.len());

            for (i, base) in ref_seq[start_idx..end_idx].chars().enumerate() {
                let x = rect.left() + (i as f32 / range_width as f32) * rect.width();
                let base_width = rect.width() / range_width as f32;

                // Background color for base
                let bg_color = base_to_color(base);
                painter.rect_filled(
                    Rect::from_min_size(
                        Pos2::new(x, rect.top()),
                        egui::vec2(base_width, row_height - 2.0),
                    ),
                    0.0,
                    bg_color,
                );

                // Draw base letter if there's room
                if base_width >= 8.0 {
                    painter.text(
                        Pos2::new(x + base_width / 2.0, rect.top() + row_height / 2.0),
                        egui::Align2::CENTER_CENTER,
                        base.to_string(),
                        egui::FontId::monospace(font_size),
                        Color32::WHITE,
                    );
                }
            }
        }

        // Draw alignment comparison (middle row) - show matches/mismatches
        let mid_y = rect.top() + row_height;
        for (i, pos) in (view_range.0..view_range.1).enumerate() {
            let x = rect.left() + (i as f32 / range_width as f32) * rect.width();
            let base_width = rect.width() / range_width as f32;

            // Find tile for this position
            let tile = tiles.iter().find(|t| pos >= t.start as u64 && pos < t.end as u64);

            let (symbol, color) = if let Some(t) = tile {
                if t.coverage == 0.0 {
                    (' ', Color32::from_gray(40))
                } else if t.identity >= 0.99 {
                    ('|', Color32::from_rgb(100, 200, 100))  // Match
                } else if t.gap_rate > 0.1 {
                    ('-', Color32::from_rgb(150, 50, 150))   // Gap/Indel
                } else {
                    ('.', Color32::from_rgb(200, 150, 50))   // Mismatch
                }
            } else {
                (' ', Color32::from_gray(40))
            };

            painter.rect_filled(
                Rect::from_min_size(
                    Pos2::new(x, mid_y),
                    egui::vec2(base_width, row_height - 2.0),
                ),
                0.0,
                color.linear_multiply(0.5),
            );

            if base_width >= 8.0 {
                painter.text(
                    Pos2::new(x + base_width / 2.0, mid_y + row_height / 2.0),
                    egui::Align2::CENTER_CENTER,
                    symbol.to_string(),
                    egui::FontId::monospace(font_size),
                    color,
                );
            }
        }

        // Draw query sequence (bottom row)
        if let Some(q_seq) = query_seq {
            let bottom_y = rect.top() + 2.0 * row_height;
            let start_idx = view_range.0 as usize;
            let end_idx = (view_range.1 as usize).min(q_seq.len());

            for (i, base) in q_seq[start_idx..end_idx].chars().enumerate() {
                let x = rect.left() + (i as f32 / range_width as f32) * rect.width();
                let base_width = rect.width() / range_width as f32;

                let bg_color = base_to_color(base);
                painter.rect_filled(
                    Rect::from_min_size(
                        Pos2::new(x, bottom_y),
                        egui::vec2(base_width, row_height - 2.0),
                    ),
                    0.0,
                    bg_color,
                );

                if base_width >= 8.0 {
                    painter.text(
                        Pos2::new(x + base_width / 2.0, bottom_y + row_height / 2.0),
                        egui::Align2::CENTER_CENTER,
                        base.to_string(),
                        egui::FontId::monospace(font_size),
                        Color32::WHITE,
                    );
                }
            }
        }

        // Row labels
        painter.text(
            Pos2::new(rect.left() + 2.0, rect.top() + 2.0),
            egui::Align2::LEFT_TOP,
            "Ref",
            egui::FontId::proportional(9.0),
            Color32::from_gray(180),
        );
        painter.text(
            Pos2::new(rect.left() + 2.0, rect.top() + row_height + 2.0),
            egui::Align2::LEFT_TOP,
            "Aln",
            egui::FontId::proportional(9.0),
            Color32::from_gray(180),
        );
        painter.text(
            Pos2::new(rect.left() + 2.0, rect.top() + 2.0 * row_height + 2.0),
            egui::Align2::LEFT_TOP,
            "Qry",
            egui::FontId::proportional(9.0),
            Color32::from_gray(180),
        );
    } else {
        // Zoomed out - show aggregated view with indel density
        painter.text(
            Pos2::new(rect.center().x, rect.center().y),
            egui::Align2::CENTER_CENTER,
            format!("Zoom in to see bases ({:.0} bp/px)", bases_per_pixel),
            egui::FontId::proportional(12.0),
            Color32::from_gray(120),
        );

        // Still show indel density as colored bars
        for tile in tiles {
            if tile.coverage == 0.0 {
                continue;
            }

            let x_start = rect.left() + (tile.start as f32 - view_range.0 as f32) / range_width as f32 * rect.width();
            let x_end = rect.left() + (tile.end as f32 - view_range.0 as f32) / range_width as f32 * rect.width();

            // Show insertions (top) and deletions (bottom) density
            let ins_height = (tile.insertion_count as f32 / 10.0).min(1.0) * (height / 3.0);
            let del_height = (tile.deletion_count as f32 / 10.0).min(1.0) * (height / 3.0);

            if ins_height > 1.0 {
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(x_start, rect.top()),
                        Pos2::new(x_end.max(x_start + 1.0), rect.top() + ins_height),
                    ),
                    0.0,
                    Color32::from_rgb(100, 150, 255),  // Blue for insertions
                );
            }

            if del_height > 1.0 {
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(x_start, rect.bottom() - del_height),
                        Pos2::new(x_end.max(x_start + 1.0), rect.bottom()),
                    ),
                    0.0,
                    Color32::from_rgb(255, 100, 100),  // Red for deletions
                );
            }
        }

        // Legend
        painter.text(
            Pos2::new(rect.right() - 60.0, rect.top() + 5.0),
            egui::Align2::LEFT_TOP,
            "Ins ▲",
            egui::FontId::proportional(9.0),
            Color32::from_rgb(100, 150, 255),
        );
        painter.text(
            Pos2::new(rect.right() - 60.0, rect.bottom() - 15.0),
            egui::Align2::LEFT_BOTTOM,
            "Del ▼",
            egui::FontId::proportional(9.0),
            Color32::from_rgb(255, 100, 100),
        );
    }
}

