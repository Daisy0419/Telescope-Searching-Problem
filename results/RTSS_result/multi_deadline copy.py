import matplotlib.pyplot as plt
import pandas as pd

# Load your CSV
file = "out_GW200216_220804.csv"
# GW191216_213338_7dt //1223
# GW191230_180458_7dt 1023
# GW200105_162426_7dt 5780
# filtered_GW200112_155838_7dt //3716
# filtered_GW200220_061928_7dt #3607
# filtered_GW200216_220804_7dt #2303
# filtered_GW200302_015811_7dt 6156
# filtered_GW200322_091133_7dt_separate 3930

df = pd.read_csv(file)

# Define beta weights for each deadline
betas = {
    50: 1.0,
    300: 0.5,
    600: 0.3,
    800: 0.1
}

# Filter GCP method only
df_gcp = df[df["Method"] == "Hoogeveen"]

# print("Available Budgets (GCP):", df[df["Method"] == "GCP"]["Budget"].unique())
# print("Beta deadlines:", list(betas.keys()))


# Create lists for matched deadlines, FoM, and discounted FoM
matched_deadlines = []
fom_values = []
discounted_scores = []

# Match deadlines to FoM and apply weighting
for D, beta in betas.items():
    row = df_gcp[df_gcp["Budget"] == D]
    if not row.empty:
        fom = row["SumProb"].values[0]
        matched_deadlines.append(D)
        fom_values.append(fom)
        discounted_scores.append(beta * fom)

# Find the best selection
best_idx = max(range(len(discounted_scores)), key=lambda i: discounted_scores[i])

# Plot
plt.figure(figsize=(6, 3), dpi=300)
# plt.plot(matched_deadlines, discounted_scores, marker='o', label='Weighted FoM (β × FoM)')
plt.plot(matched_deadlines, discounted_scores, marker='o')
plt.scatter(matched_deadlines[best_idx], discounted_scores[best_idx],
            color='red', zorder=5, label=f'Selected (D={matched_deadlines[best_idx]})')
plt.xlabel('Deadline (s)', fontsize=12)
plt.ylabel('Expected FoM', fontsize=12)
# plt.title('Deadline-Aware Prioritization')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('discounted_reward_vs_deadline.png', dpi=300)
plt.show()
