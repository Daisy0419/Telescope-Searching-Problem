import pandas as pd
import numpy as np
import plotly.graph_objects as go

def visualize_tiles(df, path_ranks=[], method_name=""):
    coords = list(zip(df['Rank'], df['RA'], df['Dec'], df['Probability']))

    ranks = []
    probs = []
    x_points = []
    y_points = []
    z_points = []

    for rank, ra_deg, dec_deg, prob in coords:
        ra_rad = np.deg2rad(ra_deg)
        dec_rad = np.deg2rad(dec_deg)
        theta = np.pi/2 - dec_rad 
        phi = ra_rad

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        ranks.append(rank)
        probs.append(prob)
        x_points.append(x)
        y_points.append(y)
        z_points.append(z)

    rank_to_xyzp = {rank: (x, y, z, prob) for rank, x, y, z, prob in zip(ranks, x_points, y_points, z_points, probs)}

    fig = go.Figure()

    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_sphere = np.outer(np.cos(u), np.sin(v))
    y_sphere = np.outer(np.sin(u), np.sin(v))
    z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))

    fig.add_trace(go.Surface(
        x=x_sphere, y=y_sphere, z=z_sphere,
        opacity=0.2, colorscale=[[0, 'gray'], [1, 'gray']],
        showscale=False
    ))

    fig.add_trace(go.Scatter3d(
        x=x_points,
        y=y_points,
        z=z_points,
        mode='markers',
        marker=dict(
            size=10,
            color=probs,
            colorscale='Viridis',
            colorbar=dict(title="Probability"),
            opacity=0.8
        ),
        hovertemplate='Rank: %{customdata[0]}<br>Prob: %{customdata[1]:.6f}',
        customdata=np.array(list(zip(ranks, probs)))
    ))

    if path_ranks:
        path_x = []
        path_y = []
        path_z = []
        for rank in path_ranks[1:-1]:
            x, y, z, _ = rank_to_xyzp[rank]
            path_x.append(x)
            path_y.append(y)
            path_z.append(z)

        fig.add_trace(go.Scatter3d(
            x=path_x,
            y=path_y,
            z=path_z,
            mode='lines+markers',
            line=dict(color='red', width=3),
            marker=dict(size=4, color='red'),
            name=f'Path ({method_name})',   
            hoverinfo='none'
        ))

    fig.update_layout(
        title=f"Visualize with Original RA/Dec Points ({method_name})", 
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='cube',
            xaxis=dict(range=[-1, 1]),
            yaxis=dict(range=[-1, 1]),
            zaxis=dict(range=[-1, 1])
        ),
        width=1000,
        height=1000
    )

    fig.show()


if __name__ == "__main__":
    # data_file = '../../data_RTSS/filtered_GW191219_163120.fits_slow_deep.csv'
    # result_file = "out_deepslow_GW191219_163120.csv"

    # data_file = '../../data_RTSS/filtered_GW191219_163120.fits_7dt.csv'
    # result_file = "out_GW191219_163120.csv"

    data_file = "../../data_RTSS_0518/large/filtered_GW200112_155838_7dt.csv"
    result_file = "large/out_filtered_GW200112_155838_7dt.csv"

    # data_file = "../../data_RTSS/filtered_GW200112_155838.fits_slow_deep.csv"
    # result_file = "out_deepslow_GW200112_155838.csv"

    data_df = pd.read_csv(data_file)
    path_df = pd.read_csv(result_file)

 
    budget_wanted   = 1000

    filtered = path_df[path_df["Budget"] == budget_wanted]
    unique_methods = filtered["Method"].unique()
    print("Identified Methods:", unique_methods)

    paths = {}  

    for method in unique_methods:
        method_df = filtered[filtered["Method"] == method]
        if len(method_df) > 0:
            paths[method] = list(map(int, method_df.iloc[0]["Path"].split()))
        else:
            print(f"No row found for method {method} with these parameters.")

    print("\nAvailable Paths:", list(paths.keys()))

    for method in list(paths.keys()):
        selected_path = paths[method]
        selected_path[1]=1
        visualize_tiles(data_df, selected_path[1:], method_name=method)


