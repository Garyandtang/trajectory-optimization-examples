%the_MAIN_block_move.m
%
% This script is the main script for block move example comparison.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
%%% Setup casadi solver 
addpath(genpath("D:\software\casadi-windows-matlabR2016a-v3.5.5"))
import casadi.*

%%% provide configuration for dynamic and problem
config = setBlockMoveConfig();

%%% compute