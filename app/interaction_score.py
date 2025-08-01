def assess_risk(interaction_energy: float) -> str:
    """Assess interaction risk based on energy."""
    if interaction_energy < -0.5:
        return "🔴 High Risk (Strong Binding)"
    elif interaction_energy < -0.2:
        return "🟠 Moderate Risk"
    elif interaction_energy < 0.2:
        return "🟢 Low Risk"
    else:
        return "🟢 No Significant Interaction"