ii = 1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/horse.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/cow.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/camel.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/venus.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/armadillo_0.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/armadillo_1.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/armadillo_4.off'); ii = ii+1;
[sig(ii).freq, sig(ii).wav] = myfunction('meshes/armadillo_10.off'); ii = ii+1;

calculateDTW(sig);