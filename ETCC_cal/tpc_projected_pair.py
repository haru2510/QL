"""散乱点近傍のstripをXZ/YZ別にフィットし、2通りの3D対応を保持する。"""
import numpy as np
from numba import njit


@njit
def projected_hits(cells, vertex, radius, velocity, offset, gem):
    """3DのX×Y候補を作らず、各投影の円内にあるstripを取り出す。"""
    hits = []
    for axis, clock, coarse, bits in cells:
        dz = velocity * (clock - offset) + gem - vertex[2]
        if dz < 3:
            continue
        for bit in range(8):
            if bits & (1 << bit):
                transverse = -153.4 + .4 * (8 * coarse + bit) - vertex[axis]
                if transverse**2 + dz**2 <= radius**2:
                    hits.append((float(axis), dz, transverse))
    return hits


@njit
def two_slopes(hits):
    """頂点を通る各枝を角度投票で初期化し、stripの割当てと直線fitを反復する。"""
    slopes = np.full((2, 2), np.nan)
    for axis in range(2):
        points = hits[hits[:, 0] == axis]
        if len(points) < 12:
            continue
        z, coordinate = points[:, 1], points[:, 2]
        angles = np.arctan2(coordinate, z)
        residual = np.empty(len(z))
        directions = np.empty(2)
        available = np.ones(len(z), np.bool_)
        for branch in range(2):
            best_count, best_angle = 0, 0.
            for grid in range(-80, 81, 2):
                angle = grid * np.pi / 180.
                distance = np.abs(coordinate*np.cos(angle) - z*np.sin(angle))
                count = np.sum(available & (distance <= 1.6))
                if count > best_count:
                    best_count, best_angle = count, angle
            directions[branch] = best_angle
            residual[:] = np.abs(coordinate*np.cos(best_angle) - z*np.sin(best_angle))
            available &= residual > 1.6
        valid = True
        for _ in range(5):
            d0 = np.abs(coordinate*np.cos(directions[0])-z*np.sin(directions[0]))
            d1 = np.abs(coordinate*np.cos(directions[1])-z*np.sin(directions[1]))
            for branch in range(2):
                mask = ((d0 <= d1) & (d0 <= 1.6) if branch == 0
                        else (d1 < d0) & (d1 <= 1.6))
                if mask.sum() < 6 or np.ptp(z[mask]) < 5.:
                    valid = False
                    break
                zz, xx = z[mask], coordinate[mask]
                directions[branch] = np.arctan(np.sum(zz*xx)/np.sum(zz**2))
            if not valid:
                break
        if valid and abs(directions[0]-directions[1]) >= 8*np.pi/180:
            slopes[axis] = np.sort(np.tan(directions))
    return slopes


def pairing_angles(slopes, x_reference, y_reference):
    """共通ΔZで枝を比較し、対応A/Bのトラック間方向をビーム横断面へ射影する。"""
    dx, dy = slopes[:, 1] - slopes[:, 0]
    vectors = np.array([[dx, dy, 0.], [dx, -dy, 0.]])
    u, v = vectors @ x_reference, vectors @ y_reference
    return np.where(u*u+v*v > 1e-16, np.degrees(np.arctan2(v,u)) % 180, np.nan)
