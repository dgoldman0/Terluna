## Training

### Round 1

Initial training was done on a 150 element sample of synthetic generated conversations about the foundations of the world, such as its history and current conditions. 

Input model: gpt-4o-mini-2024-07-18
Training data: [longer_responses.jsonl](training_data/foundation/longer_responses.jsonl)
Epochs: 3, Batch size 1, LR multiplier: 1.8, Seed: 831949690
Output model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0:CNqdXsDV

Analysis: Responses were cut very short due to overfitting on the sample and the length of each entry in each sample being short. 

### Round 2

Constructed synthetic data set for varied length assistant outputs.

Input model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0:CNqdXsDV
Training data: [sample_corrections.jsonl](training_data/foundation/correction_mode/sample_corrections.jsonl)
Epochs: 3, Batch size: 1, LR multiplier: 1.8, Seed: 1903352333
Output model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0-lengthened:CO1otnQs
