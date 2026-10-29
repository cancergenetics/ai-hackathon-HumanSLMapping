# LESSONS — process rules for running this project efficiently

Each rule: one line, with its evidence. Pruned at each milestone.

1. Check disk before copying data; read large upstream data in place behind a write guard. (Evidence: 38 GB free vs 11 GB upstream processed data, 2026-10-29.)
