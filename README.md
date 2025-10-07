## Training

### Round 1

Initial training was done on a 150 element sample of synthetic generated conversations about the foundations of the world, such as its history and current conditions. 

Input model: gpt-4o-mini-2024-07-18

Training data: [terluna_conversations_150_fixed_with_system](training_data/foundation/terluna_conversations_150_fixed_with_system.jsonl)

Epochs: 3, Batch size 1, LR multiplier: 1.8, Seed: 831949690

Output model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0:CNqdXsDV

Analysis: Responses were cut very short due to overfitting on the sample and the length of each entry in each sample being short. 

### Round 2

Constructed synthetic data set for varied length assistant outputs.

Input model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0:CNqdXsDV

Training data: [longer_responses.jsonl](training_data/foundation/longer_responses.jsonl)

Epochs: 3, Batch size: 1, LR multiplier: 1.8, Seed: 1903352333

Output model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0-lengthened:CO1otnQs

Analysis: Far from perfect, but provides much longer responses and doesn't smash unrelated points together into a meaningless paragraph as much.

### Round 3

Constructed synthetic data by going back and forth with round 2 model to identify how it responded to various inqueries, and by providing it with various corrections to align more with the desired world craft. Designed as correction mode.

Input model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0-lengthened:CO1otnQs

Training data: [correction_samples.jsonl](training_data/foundation/correction_mode/correction_samples.jsonl)

Epochs: 3, Batch size: 1, LR multiplier: 1.8, Seed: 1767361027

Output model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0-correction-mode-added:CO6Xn6Ii

Analysis: Provides some reasonable feedback and resistence against undesirable inputs. 

### Round 4

Adding more detailed history, using a combination back and forth with the round 3 model in regular and corrective mode, then transferring to agent for synthetic data generation.

Input model: ft:gpt-4o-mini-2024-07-18:personal:terluna-0-correction-mode-added:CO6Xn6Ii

Training data: [rough_draft.jsonl](training_data/history/rough_draft.jsonl)

Epochs: 3, Batch size: 1, LR multiplier: 1.8, Seed: 1014069983

Output model: 
