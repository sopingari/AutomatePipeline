import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from collections import defaultdict

def parse_pif_as_bounding_spheres(filepath):
    """
    Reads a PIFF file where each line is a single-voxel bounding box.
    Groups lines by cell_id, computing a bounding box for each cell.
    Returns two things:
      1) A dict of {cell_id: (center_x, center_y, center_z, radius)} for Body cells
      2) A single (center_x, center_y, center_z, radius) for the combined Wall vacuole
    """

    # Data structure to accumulate bounding boxes by cell_id
    # We'll store: cell_data[cell_id] = {
    #   "type": <"Body" or "Wall">,
    #   "minX": ...,
    #   "maxX": ...,
    #   "minY": ...,
    #   "maxY": ...,
    #   "minZ": ...,
    #   "maxZ": ...
    # }
    cell_data = {}

    # We'll also track an overall bounding box for all "Wall" lines
    # so we can make one big vacuole sphere
    vacuole_bounds = {
        "minX": float('inf'),
        "maxX": float('-inf'),
        "minY": float('inf'),
        "maxY": float('-inf'),
        "minZ": float('inf'),
        "maxZ": float('-inf'),
    }

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:
                # skip incomplete lines
                continue

            cell_id  = int(parts[0])
            cell_type = parts[1]  # "Body" or "Wall" (or maybe "Vacuole"?)
            x_min = int(parts[2])
            x_max = int(parts[3])
            y_min = int(parts[4])
            y_max = int(parts[5])
            z_min = int(parts[6])
            z_max = int(parts[7])

            # If this cell_id isn't in cell_data yet, initialize
            if cell_id not in cell_data:
                cell_data[cell_id] = {
                    "type": cell_type,
                    "minX": x_min,
                    "maxX": x_max,
                    "minY": y_min,
                    "maxY": y_max,
                    "minZ": z_min,
                    "maxZ": z_max
                }
            else:
                # update bounding box
                cell_data[cell_id]["minX"] = min(cell_data[cell_id]["minX"], x_min)
                cell_data[cell_id]["maxX"] = max(cell_data[cell_id]["maxX"], x_max)
                cell_data[cell_id]["minY"] = min(cell_data[cell_id]["minY"], y_min)
                cell_data[cell_id]["maxY"] = max(cell_data[cell_id]["maxY"], y_max)
                cell_data[cell_id]["minZ"] = min(cell_data[cell_id]["minZ"], z_min)
                cell_data[cell_id]["maxZ"] = max(cell_data[cell_id]["maxZ"], z_max)

            # If this line is "Wall", update the vacuole bounding box
            if cell_type == "Wall":
                vacuole_bounds["minX"] = min(vacuole_bounds["minX"], x_min)
                vacuole_bounds["maxX"] = max(vacuole_bounds["maxX"], x_max)
                vacuole_bounds["minY"] = min(vacuole_bounds["minY"], y_min)
                vacuole_bounds["maxY"] = max(vacuole_bounds["maxY"], y_max)
                vacuole_bounds["minZ"] = min(vacuole_bounds["minZ"], z_min)
                vacuole_bounds["maxZ"] = max(vacuole_bounds["maxZ"], z_max)

    # Now compute bounding spheres
    body_spheres = {}
    for cid, info in cell_data.items():
        if info["type"] == "Body":
            # compute bounding sphere
            cx = (info["minX"] + info["maxX"]) / 2.0
            cy = (info["minY"] + info["maxY"]) / 2.0
            cz = (info["minZ"] + info["maxZ"]) / 2.0
            dx = info["maxX"] - info["minX"]
            dy = info["maxY"] - info["minY"]
            dz = info["maxZ"] - info["minZ"]
            radius = 0.5 * max(dx, dy, dz)
            body_spheres[cid] = (cx, cy, cz, radius)

    # Combine all "Wall" lines into one bounding sphere
    dx = vacuole_bounds["maxX"] - vacuole_bounds["minX"]
    dy = vacuole_bounds["maxY"] - vacuole_bounds["minY"]
    dz = vacuole_bounds["maxZ"] - vacuole_bounds["minZ"]
    vac_cx = (vacuole_bounds["minX"] + vacuole_bounds["maxX"]) / 2.0
    vac_cy = (vacuole_bounds["minY"] + vacuole_bounds["maxY"]) / 2.0
    vac_cz = (vacuole_bounds["minZ"] + vacuole_bounds["maxZ"]) / 2.0
    vac_radius = 0.5 * max(dx, dy, dz)

    return body_spheres, (vac_cx, vac_cy, vac_cz, vac_radius)

def plot_spheres(body_spheres, vac_sphere, scale_factor=1.0):
    """
    Plots all body_spheres in red, and the vac_sphere in blue wireframe.
    Each sphere is scaled by 'scale_factor' to appear bigger if needed.
    """

    # Helper function to draw one sphere
    def draw_sphere(ax, center, radius, color='r', alpha=0.6, wireframe=False):
        u = np.linspace(0, 2*np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radius * np.outer(np.ones(u.size), np.cos(v))
        if wireframe:
            ax.plot_wireframe(x, y, z, color=color, alpha=alpha, linewidth=0.5)
        else:
            ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # 1) Plot vacuole sphere in blue wireframe
    vcx, vcy, vcz, vr = vac_sphere
    vcx *= scale_factor
    vcy *= scale_factor
    vcz *= scale_factor
    vr  *= scale_factor
    draw_sphere(ax, (vcx, vcy, vcz), vr, color='blue', alpha=0.2, wireframe=True)

    # 2) Plot body spheres in red
    for cid, (cx, cy, cz, r) in body_spheres.items():
        cx *= scale_factor
        cy *= scale_factor
        cz *= scale_factor
        r  *= scale_factor
        draw_sphere(ax, (cx, cy, cz), r, color='red', alpha=0.8, wireframe=False)

    # Adjust axis limits
    # Collect all coordinates
    coords = []
    coords.extend([vcx - vr, vcy - vr, vcz - vr,
                   vcx + vr, vcy + vr, vcz + vr])
    for cid, (cx, cy, cz, r) in body_spheres.items():
        coords.extend([cx*scale_factor - r*scale_factor,
                       cy*scale_factor - r*scale_factor,
                       cz*scale_factor - r*scale_factor,
                       cx*scale_factor + r*scale_factor,
                       cy*scale_factor + r*scale_factor,
                       cz*scale_factor + r*scale_factor])

    if coords:
        min_val = min(coords)
        max_val = max(coords)
    else:
        min_val = 0
        max_val = 1
    margin = 0.05 * (max_val - min_val)
    min_val -= margin
    max_val += margin

    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_zlim(min_val, max_val)
    ax.set_box_aspect((1,1,1))

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("After 500 steps")

    plt.show()

def main():
    pif_file = "./Output/10_24Simulation500.piff"  # adjust as needed

    # 1) Parse the PIF into bounding spheres
    body_spheres, vac_sphere = parse_pif_as_bounding_spheres(pif_file)

    # 2) Plot them
    # Increase scale_factor if you want them to appear bigger visually
    plot_spheres(body_spheres, vac_sphere, scale_factor=1.0)

if __name__ == "__main__":
    main()
