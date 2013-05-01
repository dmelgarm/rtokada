%degreelen.m
%This takes in the latitude and computes the length of 1 degree of latitude
%and longitude
%lat is entered in degrees

function [lenlat,lenlon] = degreelen(lat)
format long

a = 6378137;
b = 6356752.3142;
pic = pi/180;
lp = lat.*pic;

n = (a^2) ./ ((a.*cos(lp)).^2 + (b.*sin(lp)).^2).^0.5;
m = ((a*b)^2) ./ ((a.*cos(lp)).^2 + (b.*sin(lp)).^2).^1.5;

lenlon = pic.*cos(lp).*n;
lenlat = pic.*m;
