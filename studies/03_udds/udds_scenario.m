function S = udds_scenario()
%UDDS_SCENARIO Predeclared representative interval and comparison settings.
S=build_udds_scenario();
k=find(strcmp({S.segments.name},'iv208_mt13'));
assert(isscalar(k)); S.segments=S.segments(k);
S.comparison.filter_wn=[16 32];
S.comparison.battery_match_tol_As=0.5;
S.comparison.N_iter=[50 200];
S.comparison.horizons={'unif10','nonunif5'};
end
