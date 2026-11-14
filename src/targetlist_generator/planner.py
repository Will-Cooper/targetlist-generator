"""Greedy night plan: one block per target (target, calibrations, standard) in time order at low airmass."""
import math

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.time import Time

from . import sky, standards

WAIT_COST = 0.05  # wait a step rather than observe a target this far (airmass) above its best tonight
PRIORITY_STEP = 0.1  # airmass one priority level is worth when targets compete (1 = highest, 5 = lowest)


def steps(minutes: float) -> int:
    return math.ceil(minutes / sky.STEP.to_value(u.min) - 1e-9)


def plan_night(targets: pd.DataFrame, times: Time, location, max_airmass: float, target_min: np.ndarray,
               cal_min: float = 5.0, std_min: float = 10.0, stds: pd.DataFrame | None = None,
               std_tol: float = 0.1, std_max_sep: float = 20.0,
               priority: np.ndarray | None = None) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Schedule targets through the night ``times`` (a regular ``sky.STEP`` grid of dark time).
    From each free grid point, among targets whose whole exposure stays at airmass <= ``max_airmass`` and whose
    block ends before dawn, take the one closest to its best airmass tonight, each ``priority`` level (1 highest,
    5 lowest) counting as ``PRIORITY_STEP`` airmass; wait a step instead when the choice is
    still more than ``WAIT_COST`` above its best, still rising, and the night has spare time. Each target gets the
    standard chosen by ``standards.pick`` at the time it would be observed. No target is dropped: those that do not
    fit the packed night are kept as 'extra', centred on their best airmass, so the observer can choose on the night.
    Returns the plan (one row per target, status planned/extra, sorted by start time) and the paired standards
    (one row per plan row), or None without ``stds``."""
    if np.any(~np.isfinite(target_min)) or np.any(target_min <= 0):
        raise ValueError('Target times must be positive numbers of minutes')
    rank = np.zeros(len(targets)) if priority is None else (np.asarray(priority) - 1) * PRIORITY_STEP
    am = sky.airmass(targets['ra'].to_numpy(), targets['dec'].to_numpy(), times, location)  # (T, N)
    best = am.min(axis=0)
    n_exp = np.array([steps(m) for m in target_min])
    n_block = np.array([steps(m + cal_min + std_min) for m in target_min])
    n_std_start = np.array([steps(m + cal_min) for m in target_min])
    left, k, order, starts = set(range(len(targets))), 0, [], []
    while left and k < len(times):
        costs = {i: am[k:k + n_exp[i], i].mean() - best[i] for i in left
                 if k + n_block[i] <= len(times) and (am[k:k + n_exp[i], i] <= max_airmass).all()}
        if not costs:
            k += 1
            continue
        pick = min(costs, key=lambda i: costs[i] + rank[i])
        spare = len(times) - k - n_block[list(left)].sum()
        rising = am[k + 1:k + 1 + n_exp[pick], pick].mean() < am[k:k + n_exp[pick], pick].mean()
        if costs[pick] > WAIT_COST and rising and spare > 0:  # waiting improves it and the night has room
            k += 1
            continue
        order.append(pick)
        starts.append(k)
        left.remove(pick)
        k += n_block[pick]
    status = ['planned'] * len(order)
    for i in sorted(left):  # no free slot: keep it anyway, centred on its best airmass, for choosing on the night
        order.append(i)
        starts.append(int(np.clip(am[:, i].argmin() - n_exp[i] // 2, 0, len(times) - n_exp[i])))
        status.append('extra')

    idx, k0 = np.array(order, dtype=int), np.array(starts, dtype=int)
    by_time = np.argsort(k0, kind='stable')
    idx, k0, status = idx[by_time], k0[by_time], np.array(status)[by_time]
    iso = np.array([t[:16] for t in times.iso])
    mean_am = np.array([am[s:s + n_exp[i], i].mean() for i, s in zip(idx, k0, strict=True)])
    plan = targets.iloc[idx].reset_index(drop=True).assign(
        status=status, target_start=iso[k0],
        target_end=[(times[s] + target_min[i] * u.min).iso[:16] for i, s in zip(idx, k0, strict=True)],
        mean_airmass=mean_am.round(3))
    if stds is None:
        return plan, None
    std_k = np.minimum(k0 + n_std_start[idx], len(times) - 1)
    std_at = np.minimum(std_k + steps(std_min) // 2, len(times) - 1)
    paired, notes, std_am, sep = standards.pick(plan, stds, times, location, std_at, mean_am, std_tol, std_max_sep)
    return plan.assign(standard=paired['shortname'].to_numpy(), std_start=iso[std_k], std_airmass=std_am.round(3),
                       std_sep=sep.round(1), std_note=notes), paired
