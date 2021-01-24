function [ line_path ] = registerSP( AA_1_coords, AA_2_coords, SP_mask )
%UNTITLED2 Summary ofe this function goes here
%   Detailed explanation goes here
eroded = bwneighborerode(SP_mask,5,20);
trim_width = 10;
eproc = trimArray( eroded,trim_width,SP_mask );
allowable_mask = bwmorph(eproc,'thin',5);
allowable_mask = trimArray( allowable_mask,trim_width,SP_mask );
dists1 = bwdistgeodesic(allowable_mask,AA_1_coords(1),AA_1_coords(2));
dists1(isnan(dists1)) = inf;
dists2 = bwdistgeodesic(allowable_mask,AA_2_coords(1),AA_2_coords(2));
dists2(isnan(dists2)) = inf;
dists = dists1 + dists2;
dists = round(dists * 8) / 8;

paths = imregionalmin(dists);
line_path = sparse(bwmorph(paths,'thin',inf));

end

