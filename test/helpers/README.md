Shared test support lives here: optional-companion loaders, native-probe
reporting, reference-system sandboxes and scenarios reused by parent and child
processes. Helpers are included by test drivers; they are not standalone suites.

Keep a helper beside its test when only that test uses it. Move it here when
separate tests or subprocesses need the same setup or scenario.
