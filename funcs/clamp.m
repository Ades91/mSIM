function out = clamp(in,vmin,vmax)

out = in;
map = in > vmax;
out(map) = vmax;

map = in < vmin;
out(map) = vmin;