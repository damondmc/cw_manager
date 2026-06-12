import numpy as np
from astropy.io import fits

"""
Coordinate transforms between physical sky (RA/Dec) and Weave's
reduced supersky (RSSKY) parameterisation.

References
----------
LALSuite SuperSkyMetrics.c (SSM_AlignedToReduced / SSM_ReducedToAligned):
https://lscsoft.docs.ligo.org/lalsuite/7.26/lalpulsar/_supersky_metrics_8c_source.html
"""


def sky_unit_vector(alpha, delta):
    """Return the unit vector."""
    cos_delta = np.cos(delta)
    return np.array(
        [cos_delta * np.cos(alpha), cos_delta * np.sin(alpha), np.sin(delta)]
    )


def load_align_sky(setup_file):
    """Read the align_sky rotation matrix from a Weave setup FITS file."""
    with fits.open(setup_file) as hdul:
        data = hdul["semi_rssky_transf"].data["align_sky"][0]
    assert data.size == 9, f"align_sky has {data.size} elements, expected 9"
    return np.array(data).reshape(3, 3)


def get_sign_C(alpha, delta, R_asky):
    """Return sign(asky[2]) for a physical sky point — useful as a sanity check."""
    assert R_asky.shape == (3, 3), f"align_sky must be (3,3), got {R_asky.shape}"
    return float(np.sign((R_asky @ sky_unit_vector(alpha, delta))[2]))


def physical_to_rssky(alpha, delta, R_asky):
    """
    Convert physical (alpha, delta) [rad] to reduced supersky coordinates (sky_A, sky_B).
    """
    R_asky = np.asarray(R_asky).reshape(3, 3)
    n = sky_unit_vector(alpha, delta)
    na = R_asky @ n  # aligned sky coordinates
    r = np.linalg.norm(na, axis=0)
    rss_A = np.sign(na[2]) * (na[0] / r + 1.0)
    rss_B = na[1] / r
    return rss_A, rss_B


def rssky_to_physical(rss_A, rss_B, R_asky):
    """
    Convert reduced supersky coordinates (rss_A, rss_B) to physical (alpha, delta) [rad].

    Note: sign(rss_A) encodes the hemisphere — no separate sign_C needed.

    Parameters
    ----------
    rss_A, rss_B : float or ndarray
    R_asky : ndarray, shape (9,) or (3,3)
        Rotation matrix from Weave setup file semi_rssky_transf for aligned sky coordinates.
    """
    R_asky = np.asarray(R_asky).reshape(3, 3)
    rss_A = np.asarray(rss_A, dtype=float)
    rss_B = np.asarray(rss_B, dtype=float)
    dA = np.fabs(rss_A) - 1.0
    Rmax = np.maximum(1.0, np.sqrt(dA**2 + rss_B**2))
    n_A = dA / Rmax
    n_B = rss_B / Rmax
    n_C = np.sign(rss_A) * np.sqrt(np.maximum(0.0, 1.0 - n_A**2 - n_B**2))
    na = np.stack([n_A, n_B, n_C], axis=-1)  # shape (..., 3) #
    # ssky = R_asky^T * asky  (i.e. n = R^T @ na for column vectors)
    n_phys = na @ R_asky
    rho = np.sqrt(n_phys[..., 0] ** 2 + n_phys[..., 1] ** 2)
    delta = np.arctan2(n_phys[..., 2], rho)
    alpha = np.arctan2(n_phys[..., 1], n_phys[..., 0]) % (2 * np.pi)
    return alpha, delta


def rssky_offsets_to_physical(rss_A_0, rss_B_0, d_rss_A, d_rss_B, R_asky):
    """
    Convert RSSKY grid offsets to physical (d_alpha, d_delta) offsets.

    Uses the full exact inverse transform (not the linearized Jacobian).

    Parameters
    ----------
    rss_A_0, rss_B_0 : float
        RSSKY coordinates of the reference (target) sky point.
    d_rss_A, d_rss_B : float or ndarray
        RSSKY offsets (from sky_hex_offsets using RSSKY spacings).
    R_asky : ndarray
        Rotation matrix from Weave setup file for aligned sky coordinates.

    Returns
    -------
    d_alpha, d_delta : ndarray
        Physical RA/Dec offsets in radians.
    """
    alpha_0, delta_0 = rssky_to_physical(rss_A_0, rss_B_0, R_asky)
    alpha, delta = rssky_to_physical(
        rss_A_0 + np.asarray(d_rss_A),
        rss_B_0 + np.asarray(d_rss_B),
        R_asky,
    )
    return alpha - alpha_0, delta - delta_0


def read_sky_grid_info(weave_result_file):
    """
    Read sky grid info from a raw Weave result FITS primary header.

    Keys used:
        NSEMITMPL SSKYA  — raw template count along Weave's sky_A axis
                           (= this module's rss_B direction)
        NSEMITMPL SSKYB  — *cumulative* template count (count_A * count_B);
                           dividing by NSEMITMPL SSKYA gives the raw count
                           along Weave's sky_B axis (= this module's rss_A direction)
        PROGARG ALPHA    — physical RA range "min,max" [rad]
        PROGARG DELTA    — physical Dec range "min,max" [rad]

    Returns
    -------
    n_rssA : int
    n_rssB : int
    alpha_range : (float, float)   (min, max) RA  [rad]
    delta_range : (float, float)   (min, max) Dec [rad]
    """
    with fits.open(weave_result_file) as hdul:
        h = hdul[0].header
    n_rssB = int(h["NSEMITMPL SSKYA"])
    n_rssA = int(h["NSEMITMPL SSKYB"]) // n_rssB
    alpha_min, alpha_max = (float(x) for x in h["PROGARG ALPHA"].split(","))
    delta_min, delta_max = (float(x) for x in h["PROGARG DELTA"].split(","))
    return n_rssA, n_rssB, (alpha_min, alpha_max), (delta_min, delta_max)


def generate_sky_grid(alpha_range, delta_range, n_rssA, n_rssB, R_asky):
    """
    Generate physical (alpha, delta) grid points on the hexagonal (triangular /
    A_2*-type) covering lattice that lalpulsar_Weave's LatticeTiling places in
    RSSKY space.

    Within each row (fixed rss_B) templates are evenly spaced along rss_A;
    consecutive rows are laterally shifted by half the within-row spacing.
    This row-to-row shear is exact for a hexagonal lattice whenever the
    sky-sky metric block is diagonal in (rss_A, rss_B) -- true here, since
    semi_rssky_metric has M_AB = 0 -- and is what makes the real grid look
    "sheared/rotated" relative to a plain axis-aligned rectangle (verified
    against an actual Weave search grid: mean nearest-neighbour residual
    drops from ~1.4e-5 rad for a plain meshgrid to ~5e-6 rad with the shear).

    The RSSKY extents (hence row/column spacings) are estimated by projecting
    the physical RA/Dec range edges, same as before; the shear is centered on
    the middle row so the grid stays anchored at the range center.

    Parameters
    ----------
    alpha_range : (float, float)   RA range  (min, max) [rad]
    delta_range : (float, float)   Dec range (min, max) [rad]
    n_rssA, n_rssB : int                 Template counts along rss_A and rss_B
    R_asky : ndarray, shape (3, 3) Rotation matrix from load_align_sky()

    Returns
    -------
    alphas, deltas : ndarray, shape (n_rssA * n_rssB,)
        Physical sky coordinates of all grid points.
    """
    assert R_asky.shape == (3, 3), f"R_asky must be (3,3), got {R_asky.shape}"

    alpha_c = 0.5 * (alpha_range[0] + alpha_range[1])
    delta_c = 0.5 * (delta_range[0] + delta_range[1])
    rss_A_0, rss_B_0 = physical_to_rssky(alpha_c, delta_c, R_asky)

    # Estimate rss_A extent (and within-row spacing) from RA edges at fixed Dec center
    if n_rssA > 1:
        rss_A_lo = physical_to_rssky(alpha_range[0], delta_c, R_asky)[0]
        rss_A_hi = physical_to_rssky(alpha_range[1], delta_c, R_asky)[0]
        rss_A_vals = np.linspace(rss_A_lo, rss_A_hi, n_rssA)
        spacing_A = (rss_A_hi - rss_A_lo) / (n_rssA - 1)
    else:
        rss_A_vals = np.array([rss_A_0])
        spacing_A = 0.0

    # Estimate rss_B extent from Dec edges at fixed RA center
    if n_rssB > 1:
        rss_B_lo = physical_to_rssky(alpha_c, delta_range[0], R_asky)[1]
        rss_B_hi = physical_to_rssky(alpha_c, delta_range[1], R_asky)[1]
        rss_B_vals = np.linspace(rss_B_lo, rss_B_hi, n_rssB)
    else:
        rss_B_vals = np.array([rss_B_0])

    # Hexagonal-lattice row shear: each row is offset from its neighbour by
    # half the within-row spacing (sign fixed empirically from a real Weave
    # grid; the 1/2 magnitude is exact for any diagonal sky-sky metric block).
    # Centering on the middle row keeps the grid anchored at the range center.
    row_offset = np.arange(n_rssB) - 0.5 * (n_rssB - 1)
    rss_A_grid = rss_A_vals[np.newaxis, :] - 0.5 * spacing_A * row_offset[:, np.newaxis]
    rss_B_grid = np.broadcast_to(rss_B_vals[:, np.newaxis], rss_A_grid.shape)

    alphas, deltas = rssky_to_physical(rss_A_grid.ravel(), rss_B_grid.ravel(), R_asky)
    return alphas, deltas


def generate_sky_grid_local(alpha_range, delta_range, n_rows, n_total):
    """
    Generate a physical (alpha, delta) sky template grid matching the row
    structure lalpulsar_Weave actually uses for *directed* searches: rows of
    (nearly) constant declination, evenly spaced across delta_range, each
    holding ~n_total / n_rows points evenly spaced across alpha_range, with
    consecutive rows laterally offset by half the within-row spacing
    (hexagonal shear).

    Unlike generate_sky_grid, this needs no align_sky / RSSKY machinery: for
    a directed search the requested sky region is small enough that physical
    (RA, Dec) IS effectively the lattice's natural frame -- confirmed against
    real Weave search grids at six widely separated sky locations, where the
    dominant nearest-neighbour direction matched the "constant declination"
    direction to within ~1-4 deg. (Assuming rows of constant RSSKY rss_B
    instead, as generate_sky_grid does, is only right by coincidence wherever
    the all-sky setup file's align_sky happens to line up physical RA/Dec
    with its rss_A/rss_B axes -- elsewhere it's off by tens of degrees.)

    n_rows and n_total should come straight from the *cumulative* counts
    NSEMITMPL SSKYA / SSKYB, NOT the average row length n_rssA = SSKYB //
    SSKYA that read_sky_grid_info derives: real row lengths vary, since the
    requested rectangle maps to a slightly tilted region in the lattice
    frame that tapers near its corners. Assuming a uniform n_total // n_rows
    per row is therefore still an approximation -- expect any residual
    mismatch to be concentrated at the grid's edges.

    Parameters
    ----------
    alpha_range : (float, float)   RA range  (min, max) [rad]
    delta_range : (float, float)   Dec range (min, max) [rad]
    n_rows : int    Number of declination rows (= NSEMITMPL SSKYA)
    n_total : int   Total number of sky templates (= NSEMITMPL SSKYB)

    Returns
    -------
    alphas, deltas : ndarray, shape (n_rows * round(n_total / n_rows),)
        Physical sky coordinates of all grid points.
    """
    delta_vals = (
        np.linspace(delta_range[0], delta_range[1], n_rows)
        if n_rows > 1
        else np.array([0.5 * (delta_range[0] + delta_range[1])])
    )

    n_per_row = max(1, round(n_total / n_rows))
    if n_per_row > 1:
        alpha_vals = np.linspace(alpha_range[0], alpha_range[1], n_per_row)
        spacing = (alpha_range[1] - alpha_range[0]) / (n_per_row - 1)
    else:
        alpha_vals = np.array([0.5 * (alpha_range[0] + alpha_range[1])])
        spacing = 0.0

    row_offset = np.arange(n_rows) - 0.5 * (n_rows - 1)
    alpha_grid = alpha_vals[np.newaxis, :] - 0.5 * spacing * row_offset[:, np.newaxis]
    delta_grid = np.broadcast_to(delta_vals[:, np.newaxis], alpha_grid.shape)

    return alpha_grid.ravel(), delta_grid.ravel()


def rssky_jacobian(rss_A_0, rss_B_0, R_asky):
    """
    Linearized Jacobian J such that [d_alpha, d_delta] = J @ [d_rss_A, d_rss_B].

    Derivation
    ----------
    rss_A = sign(asky_C) * (asky_A + 1)  →  d_asky_A = sign(rss_A) * d_rss_A

    Constrained tangent vectors in the aligned frame (unit-sphere constraint
    asky_A² + asky_B² + asky_C² = 1 fixes the third component):
        v_A = sign(rss_A) * [1,  0,  -asky_A / asky_C]
        v_B =               [0,  1,  -asky_B / asky_C]

    Map to physical frame:  t_A = R^T @ v_A,  t_B = R^T @ v_B

    Physical offsets at reference (n_x, n_y, n_z), ρ² = n_x² + n_y²:
        J[0, 0] = (n_x * t_A[1] - n_y * t_A[0]) / ρ²    (d_alpha / d_rss_A)
        J[0, 1] = (n_x * t_B[1] - n_y * t_B[0]) / ρ²    (d_alpha / d_rss_B)
        J[1, 0] = t_A[2] / cos(delta)                     (d_delta / d_rss_A)
        J[1, 1] = t_B[2] / cos(delta)                     (d_delta / d_rss_B)

    To generate sky grid via Jacobian (alternative to edge conversion):
        d_rss = inv(J) @ [d_alpha_total, d_delta_total]
        rss_A_vals = rss_A_0 + linspace(-0.5, 0.5, n_A) * d_rss[0]
        rss_B_vals = rss_B_0 + linspace(-0.5, 0.5, n_B) * d_rss[1]
    """
    R_asky = np.asarray(R_asky).reshape(3, 3)
    sign_A = float(np.sign(rss_A_0))

    # Recover aligned coordinates from RSSKY (Rmax=1 for interior points)
    n_A_0 = np.fabs(rss_A_0) - 1.0
    n_B_0 = rss_B_0
    n_C_0 = sign_A * np.sqrt(max(0.0, 1.0 - n_A_0**2 - n_B_0**2))

    # Constrained tangent vectors in the aligned frame
    v_A = sign_A * np.array([1.0, 0.0, -n_A_0 / n_C_0])
    v_B = np.array([0.0, 1.0, -n_B_0 / n_C_0])

    t_A = R_asky.T @ v_A
    t_B = R_asky.T @ v_B

    alpha_0, delta_0 = rssky_to_physical(rss_A_0, rss_B_0, R_asky)
    nx = np.cos(delta_0) * np.cos(alpha_0)
    ny = np.cos(delta_0) * np.sin(alpha_0)
    rho2 = nx**2 + ny**2

    return np.array(
        [
            [(nx * t_A[1] - ny * t_A[0]) / rho2, (nx * t_B[1] - ny * t_B[0]) / rho2],
            [t_A[2] / np.cos(delta_0), t_B[2] / np.cos(delta_0)],
        ]
    )
