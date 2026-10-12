module RuntimeServices

# Despite living in separate files (`libGRAM.dylib`'s statically-linked CSPICE
# copy vs. SPICE.jl's own libcspice), the two are not isolated at the OS
# symbol-resolution level: `nm -gU libGRAM.dylib` shows it globally exports
# CSPICE internals (`chkin_`, `chkout_`, `trcpkg_`, `subslr_c`, ...) under the
# same symbol names SPICE.jl's libcspice also exports. Concurrently invoking
# a native-GRAM call (previously serialized only by GRAM_LOCK) on one thread
# and a SpaceAGORA ephemerides/frame-transform call (serialized only by
# SPICE_LOCK) on another corrupts CSPICE's internal call-trace stack
# (SPICE(NAMESDONOTMATCH)/CHKOUT errors), because both locks were guarding
# what is actually one shared C-level critical section. GRAM_LOCK and
# SPICE_LOCK must therefore be the same lock, not merely same-shaped.
const SPICE_LOCK = ReentrantLock()
const GRAM_LOCK = SPICE_LOCK

end # module RuntimeServices
