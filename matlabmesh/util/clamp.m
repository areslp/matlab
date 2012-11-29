function [ clamped ] = clamp( val, minv, maxv )
%CLAMP Summary of this function goes here
%   Detailed explanation goes here
clamped = max(minv, min(maxv, val));
