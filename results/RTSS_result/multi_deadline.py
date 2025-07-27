import matplotlib.pyplot as plt
import pandas as pd

# Load your CSV
# path = "large"
# filtered_GW191216_213338_7dt //1223
# filtered_GW191230_180458_7dt 1023
# filtered_GW200105_162426_7dt 5780
# filtered_GW200112_155838_7dt //3716
# filtered_GW200220_061928_7dt #3607
# filtered_GW200216_220804_7dt #2303
# filtered_GW200302_015811_7dt 6156
# filtered_GW200322_091133_7dt_separate 3930
# filtered_GW200112_155838_7dt

# path = "multi_deadline"


import pandas as pd
import matplotlib.pyplot as plt

files = ["GW200112_155838", "GW200216_220804", "GW200105_162426"]

# Define beta weights for each deadline
betas = {
    50: 1.0,
    100: 0.8,
    200: 0.5,
    # 500: 0.3,
    900: 0.2,
    1200: 0.0
}
# betas = {
#     50: 1.0,
#     300: 0.5,
#     600: 0.3,
#     800: 0.1,
# }

path = "large"
plt.figure(figsize=(6, 3), dpi=300)

for fname in files:
    filepath = f"{path}/multi_filtered_{fname}_7dt.csv"  # or adjust extension path as needed
    df = pd.read_csv(filepath)
    df_gcp = df[df["Method"] == "Hoogeveen"]

    matched_deadlines = []
    discounted_scores = []
    
    for D, beta in betas.items():
        row = df_gcp[df_gcp["Budget"] == D]
        if not row.empty:
            fom = row["SumProb"].values[0]
            matched_deadlines.append(D)
            discounted_scores.append(beta * fom)

    if matched_deadlines:
        best_idx = max(range(len(discounted_scores)), key=lambda i: discounted_scores[i])
        label = fname.replace("_", "\n")  # optionally wrap text for clarity

        # Plot line
        plt.plot(matched_deadlines, discounted_scores, marker='o', label=label)
        # Highlight the best choice
        plt.scatter(matched_deadlines[best_idx], discounted_scores[best_idx],
                    color='red', zorder=5)

plt.xlabel('Deadline (s)', fontsize=12)
plt.ylabel('Expected FoM', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(title="Skymap", fontsize=9)
plt.tight_layout()
plt.savefig('multi_discounted_reward_vs_deadline.png', dpi=300)
plt.show()


# files = ["GW191216_213338", "GW200220_061928", "GW200302_015811"]

# df = pd.read_csv(file)

# # Define beta weights for each deadline
# betas = {
#     50: 1.0,
#     300: 0.5,
#     600: 0.3,
#     800: 0.1
# }

# df_gcp = df[df["Method"] == "Hoogeveen"]

# matched_deadlines = []
# fom_values = []
# discounted_scores = []

# for D, beta in betas.items():
#     row = df_gcp[df_gcp["Budget"] == D]
#     if not row.empty:
#         fom = row["SumProb"].values[0]
#         matched_deadlines.append(D)
#         fom_values.append(fom)
#         discounted_scores.append(beta * fom)

# best_idx = max(range(len(discounted_scores)), key=lambda i: discounted_scores[i])

# # Plot
# plt.figure(figsize=(6, 3), dpi=300)
# # plt.plot(matched_deadlines, discounted_scores, marker='o', label='Weighted FoM (β × FoM)')
# plt.plot(matched_deadlines, discounted_scores, marker='o')
# plt.scatter(matched_deadlines[best_idx], discounted_scores[best_idx],
#             color='red', zorder=5, label=f'Selected (D={matched_deadlines[best_idx]})')
# plt.xlabel('Deadline (s)', fontsize=12)
# plt.ylabel('Expected FoM', fontsize=12)
# # plt.title('Deadline-Aware Prioritization')
# plt.grid(True, linestyle='--', alpha=0.6)
# plt.legend(fontsize=12)
# plt.tight_layout()
# plt.savefig('discounted_reward_vs_deadline.png', dpi=300)
# plt.show()
