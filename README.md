# TrueSkill.jl

TrueSkill, the video game rating system

Julia translation of the Python [TrueSkill](https://trueskill.org) project.

## Caution

This TrueSkill project is opened under the BSD license but the
TrueSkill(TM) brand is not. Microsoft permits only Xbox Live games or
non-commercial projects to use TrueSkill(TM). If your project is
commercial, you should find another rating system.

## What is TrueSkill?

TrueSkill is a rating system among game players. It was developed by Microsoft Research and has been used on Xbox LIVE for ranking and matchmaking service. This system quantifies players’ TRUE skill points by the Bayesian inference algorithm. It also works well with any type of match rule including N:N team game or free-for-all.

This project is a Julia package which implements the TrueSkill rating system.

## Installation

```julia
using Pkg
pkg"add https://github.com/kose-y/TrueSkill.jl"
```

## Basics

### Rating, the model for skill
In TrueSkill, rating is a Gaussian distribution which starts from N(25,2532). μ is an average skill of player, and σ is a confidence of the guessed rating. A real skill of player is between μ±2σ with 95% confidence.



```julia
using TrueSkill
Rating()
```




    Rating(mu=25.0, sigma=8.333333333333334)



If some player’s rating is higher β than another player’s, the player may have about a 76% chance to beat the other player. The default value of β is 25/6.

Ratings will approach real skills through few times of the TrueSkill's Bayesian inference algorithm. How many matches TrueSkill needs to estimate real skills? It depends on the game rule. See the below table:

| Rule | Matches |
|:---|---:|
| 16P free-for-all | 3 |
| 8P free-for-all | 3 |
| 4P free-for-all | 5 |
| 2P free-for-all | 12 |
| 2:2:2:2 | 10 |
| 4:4:4:4 | 20 |
| 4:4 | 46 |
| 8:8 | 91 |


### Head-to-head (1 vs. 1) match rule


```julia
r1 = Rating()
r2 = Rating()
new_r1, new_r2 = rate(r1, r2)
```




    (Rating(mu=29.39583169299151, sigma=7.17147580700922), Rating(mu=20.604168307008482, sigma=7.17147580700922))



User should put the winner first. After the game, TrueSkill recalculates their ratings by the game rseult. You can also handle a tie game:


```julia
r1 = Rating()
r2 = Rating()
new_r1, new_r2 = rate(r1, r2; drawn=true)
```




    (Rating(mu=24.999999999999993, sigma=6.457515683245051), Rating(mu=24.999999999999993, sigma=6.457515683245051))



### Other match rules

4-player free-for-all:


```julia
r1 = Rating()
r2 = Rating()
r3 = Rating()
r4 = Rating()
new_r1, new_r2, new_r3, new_r4 = rate([r1, r2, r3, r4])
```




    4-element Array{Rating,1}:
     Rating(mu=33.20668089498382, sigma=6.348109386291753)
     Rating(mu=27.40145515734498, sigma=5.7871628096649355)
     Rating(mu=22.598544842686724, sigma=5.787162809661523)
     Rating(mu=16.793319105008695, sigma=6.348109386298544)



4-player free-for-all with ties:


```julia
r1 = Rating()
r2 = Rating()
r3 = Rating()
r4 = Rating()
new_r1, new_r2, new_r3, new_r4 = rate([r1, r2, r3, r4]; ranks=[1, 2, 2, 4])
```




    4-element Array{Rating,1}:
     Rating(mu=31.56397232349745, sigma=6.4047036446727095)
     Rating(mu=24.993092934523535, sigma=5.559362333660422)
     Rating(mu=25.006907065417195, sigma=5.559362333660494)
     Rating(mu=18.436027676561814, sigma=6.404703644728995)



2:2:2 match:


```julia
r1 = Rating()
r2 = Rating()
r3 = Rating()
r4 = Rating()
r5 = Rating()
r6 = Rating()
(r1, r2), (r3, r4), (r5, r6) = rate([[r1, r2], [r3, r4], [r5, r6]])
```




    3-element Array{Array{Rating,1},1}:
     [Rating(mu=29.72018660358779, sigma=7.541668779519052), Rating(mu=29.720186603587795, sigma=7.541668779519052)]
     [Rating(mu=25.000000000002768, sigma=7.34810769389046), Rating(mu=25.000000000002764, sigma=7.34810769389046)]
     [Rating(mu=20.279813396409434, sigma=7.541668779519555), Rating(mu=20.279813396409434, sigma=7.541668779519555)]



2:1 match:


```julia
r1 = Rating()
r2 = Rating()
r3 = Rating()
t1 = [r1]
t2 = [r2, r3]
t1, t2 = rate([t1, t2]) # t1 defeats t2
```




    2-element Array{Array{Rating,1},1}:
     [Rating(mu=33.73067114899642, sigma=7.317365362867211)]
     [Rating(mu=16.269328851003575, sigma=7.317365362867211), Rating(mu=16.269328851003575, sigma=7.317365362867211)]


