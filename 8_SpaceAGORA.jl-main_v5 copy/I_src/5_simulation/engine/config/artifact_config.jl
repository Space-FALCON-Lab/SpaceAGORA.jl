"""
    ArtifactConfig

Typed runtime configuration for result bundles, checkpoint-related output, and
warning policy around deprecated artifact behavior.
"""
Base.@kwdef struct ArtifactConfig
    save_bundle::Bool = true
    warn_deprecated_config::Bool = true
end
