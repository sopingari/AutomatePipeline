import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_vacuole_spheres(df, show_inner=False, save_path=None, scale_factor=10):
    """
    Plot spheres (APBs) inside a vacuole using data from a DataFrame
    with columns at least:
      - 'bodyType': "APB" or "Vacuole" (or "Wall")
      - 'rInner', 'rOuter' (for the vacuole row, if available)
      - 'x', 'y', 'z', 'r' (for each APB and the vacuole)
    
    The coordinates and radii are multiplied by scale_factor for visualization.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the simulation results.
    show_inner : bool, optional
        If True, also draws the vacuole's inner boundary as a wireframe.
    save_path : str, optional
        If provided, saves the figure to this file path. Otherwise, displays it.
    scale_factor : float, optional
        Factor by which to multiply the coordinates and radii for display.
    """
    
    def draw_sphere(ax, center, radius, color='blue', alpha=0.2, wireframe=False):
        u = np.linspace(0, 2*np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        # Scale the sphere coordinates
        x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radius * np.outer(np.ones(u.size), np.cos(v))
        if wireframe:
            ax.plot_wireframe(x, y, z, color=color, alpha=alpha, linewidth=0.5)
        else:
            ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Look for a row with bodyType either Vacuole or Wall
    vacuole_rows = df.loc[df['bodyType'].isin(['Vacuole', 'Wall'])]
    if vacuole_rows.empty:
        print("Warning: No vacuole (or wall) row found. Plotting APBs only.")
        vacuole = None
    else:
        vacuole = vacuole_rows.iloc[0].copy()
        # If rOuter or rInner aren't provided, use r
        if 'rOuter' not in vacuole:
            vacuole['rOuter'] = vacuole['r']
        if 'rInner' not in vacuole:
            vacuole['rInner'] = vacuole['r']
    
    # Multiply coordinates and radii by the scale factor for display
    if vacuole is not None:
        vacuole['x'] *= scale_factor
        vacuole['y'] *= scale_factor
        vacuole['z'] *= scale_factor
        vacuole['rOuter'] *= scale_factor
        vacuole['rInner'] *= scale_factor

    # Process APBs: scale coordinates and radii
    apbs = df.loc[df['bodyType'] == 'APB'].copy()
    apbs['x'] *= scale_factor
    apbs['y'] *= scale_factor
    apbs['z'] *= scale_factor
    apbs['r'] *= scale_factor

    if vacuole is not None:
        draw_sphere(ax,
                    center=(vacuole['x'], vacuole['y'], vacuole['z']),
                    radius=vacuole['rOuter'],
                    color='skyblue',
                    alpha=0.2,
                    wireframe=False)
    
        if show_inner and 'rInner' in vacuole:
            draw_sphere(ax,
                        center=(vacuole['x'], vacuole['y'], vacuole['z']),
                        radius=vacuole['rInner'],
                        color='blue',
                        alpha=0.3,
                        wireframe=True)
    
    for _, row in apbs.iterrows():
        draw_sphere(ax,
                    center=(row['x'], row['y'], row['z']),
                    radius=row['r'],
                    color='red',
                    alpha=0.8,
                    wireframe=False)
    
    # Determine axis limits based on scaled data
    coords = []
    if not apbs.empty:
        coords.extend([apbs['x'].min(), apbs['y'].min(), apbs['z'].min(),
                       apbs['x'].max(), apbs['y'].max(), apbs['z'].max()])
    if vacuole is not None:
        coords.extend([
            vacuole['x'] - vacuole['rOuter'], vacuole['x'] + vacuole['rOuter'],
            vacuole['y'] - vacuole['rOuter'], vacuole['y'] + vacuole['rOuter'],
            vacuole['z'] - vacuole['rOuter'], vacuole['z'] + vacuole['rOuter']
        ])
    if coords:
        min_val, max_val = np.nanmin(coords), np.nanmax(coords)
    else:
        min_val, max_val = 0, 1
    margin = 0.05 * (max_val - min_val)
    min_val -= margin
    max_val += margin
    
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_zlim(min_val, max_val)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Result of vacuole_gen.py')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
        plt.close(fig)
    else:
        plt.show()
