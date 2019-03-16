@echo off
taskkill /fi "WindowTitle eq Gnuplot*"
start grapher.plt