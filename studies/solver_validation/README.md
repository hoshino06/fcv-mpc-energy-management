# Solver validation

`check_solver_snapshot.m` compares DDP and fmincon from the same state, demand
preview, nominal input, initial guess, and battery constraint. These checks
support interpretation of suboptimal DDP solutions; they are separate from the
closed-loop timing study.
