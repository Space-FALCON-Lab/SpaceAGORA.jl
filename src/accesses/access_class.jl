"""
    Access

Represents a time window during which a spacecraft is accessible.

# Fields
- `start_time::Float64`: Start time of the access window (in seconds)
- `end_time::Float64`: End time of the access window (in seconds)
- `spacecraft_id::String`: Unique identifier for the spacecraft

# Constructor
    Sat_Access(start_time::Float64, end_time::Float64, target_id::Int64, spacecraft_id::Int64)
"""

struct Access
    start_time::Float64
    end_time::Float64
    target_id::Int64
    spacecraft_id::Int64

end

function start_access(start_time::Float64, target_id::Int64, spacecraft_id::Int64)
    return Access(start_time, -1.0, target_id, spacecraft_id)
end
function end_access(access::Access, end_time::Float64)
    return Access(access.start_time, end_time, access.target_id, access.spacecraft_id)
end

