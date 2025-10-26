def pick(poses, cfg: dict):
    # MVP: choose first pose path. Later: cluster and score-based selection.
    return poses[0] if poses else None
