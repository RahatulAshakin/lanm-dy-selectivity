# ProteinMPNN Round-2 Smoke

Deterministic Phase 7B ProteinMPNN smoke generation across the narrow round-2 redesign campaigns.

- ProteinMPNN root: `/home/ashak/apps/ProteinMPNN`
- num_seq_per_target: `20`
- sampling_temp: `0.1 0.15`
- seed: `37`
- batch_size: `1`

| seed_scaffold | redesignable_positions | unique_sequences | mutation_count_distribution |
| --- | --- | ---: | --- |
| am1_mex_ss_only_u02 | 18,24,27,38,39,40,44,65,69,73,75,76,79,80,83,97,98,100 | 29/40 | 3:7,4:6,5:14,6:11,7:2 |
| hans_pocket_ss_only_u02 | 18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100 | 33/40 | 3:1,4:14,5:16,6:8,7:1 |
| hans_interface_ss_plus_if_u04 | 18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100 | 12/40 | 4:14,5:23,6:3 |

## am1_mex_ss_only_u02

- Seed scaffold: `am1_mex_ss_only_u02` from `results/design_inputs/proteinmpnn_round2/campaigns/am1_mex_ss_only_u02/am1_mex_ss_only_u02.pdb`
- Designed-chain seed sequence: `VDIAAFDPDDDGTIDLKEALAAGSAAFDKLDPDKDGTLDAKELKGRVSEADLKKYDPDKDGTLDKKEYLAAVEAQFKAANPDNDGTIDAAELASPAGSALVNLIR`
- Redesignable positions: residue ids `A38,A44,A47,A56,A57,A58,A62,A83,A87,A91,A93,A94,A97,A98,A101,A115,A116,A118`; canonical `18,24,27,38,39,40,44,65,69,73,75,76,79,80,83,97,98,100`; AM1 `16,22,25,34,35,36,40,61,65,69,71,72,75,76,79,93,94,96`
- Unique sequences: `29` of `40` generated
- Mutation-count distribution vs seed scaffold: `3:7,4:6,5:14,6:11,7:2`

## hans_pocket_ss_only_u02

- Seed scaffold: `hans_pocket_ss_only_u02` from `results/design_inputs/proteinmpnn_round2/campaigns/hans_pocket_ss_only_u02/hans_pocket_ss_only_u02.pdb`
- Designed-chain seed sequence: `ASGADALKALNTDNDDSLEIAEVIHAGATTFTAINPDGDTTLESGETKGRLTEKDWARANKDGDQTLEMDEWLKILKTRFKRADANGDGKLDAAELDSKAGQGVLVMIMK`
- Redesignable positions: residue ids `A35,A41,A44,A55,A56,A57,A61,A90,A92,A93,A96,A97,A100,A110,A114,A115,A117`; canonical `18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100`; AM1 `16,22,25,34,35,36,40,69,71,72,75,76,79,89,93,94,96`
- Unique sequences: `33` of `40` generated
- Mutation-count distribution vs seed scaffold: `3:1,4:14,5:16,6:8,7:1`

## hans_interface_ss_plus_if_u04

- Seed scaffold: `hans_interface_ss_plus_if_u04` from `results/design_inputs/proteinmpnn_round2/campaigns/hans_interface_ss_plus_if_u04/hans_interface_ss_plus_if_u04.pdb`
- Designed-chain seed sequence: `ASGADALKALNSDNDDSLEIAEVIHAGATTFWAINPDGDTTLESGETKGRLTEKDWARANKDGDQTLEMDEWLKILRTRFKRADANGDGKLEAAELDSKAGQGVLVMIMK`
- Redesignable positions: residue ids `A35,A41,A44,A55,A56,A57,A61,A90,A92,A93,A96,A97,A100,A110,A114,A115,A117`; canonical `18,24,27,38,39,40,44,73,75,76,79,80,83,93,97,98,100`; AM1 `16,22,25,34,35,36,40,69,71,72,75,76,79,89,93,94,96`
- Unique sequences: `12` of `40` generated
- Mutation-count distribution vs seed scaffold: `4:14,5:23,6:3`
