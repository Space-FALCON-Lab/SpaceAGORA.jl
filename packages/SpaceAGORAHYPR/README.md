# SpaceAGORAHYPR

Optional HYPR execution for the matching SpaceAGORA checkout. Install both in a
separate Julia project with `Pkg.develop`, then load `SpaceAGORAHYPR` explicitly.
SpaceAGORA supplies the shared planner interface and public configuration types;
this companion supplies configured HYPR and robot-arm search. Both import orders
are supported.

See `docs/src/user/rpo_planner_pilot.md` in the repository for installation and
validated scope. The companion introduces no new scientific acceptance limits.
