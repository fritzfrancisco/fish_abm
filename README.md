## fish abm
Agent based model for simulating movement of fish schools.

<img src="https://github.com/fritzfrancisco/fish_abm/blob/master/screenshot_abm.png" width="500">

Image 1: First visualization of the agent based model without environment

<img src="https://github.com/fritzfrancisco/fish_abm/blob/master/screenshot_abm_environment.png" width="500">

Image 2: First visualization of the agent based model with environment


### `int main():`

- create `Mat environment` & `Mat fish`
- `create_environment()`
- `individuals fish`
- `initialize()`

  #### `void sample():`

  - counting individuals `in_zone(r = 15)`,`in_zone(r = 100)` and `in_zone(r = 200)`

    1. when alone lag is decreased [-1] and individual moves to group
    2. if `sample_rate > distrib_sample(rd)` (random number) lag is set to [5] and individual feeds. This implements that individuals by chance will sample more frequently in a "good" food patch
    3. else, if lag is 0, `sample_rate` increases [+1]
    4. else, lag is reduced [-1] and `sample_rate` increases [+1]

  #### `void feed():`

  - `sample_rate = 10 * q` (q being the quality in the quality box encountered)
  - reduce quality of quality box encountered (feeding)

  #### `void move():`

  - the direction angle is corrected back between 0 - 360 degrees and the individual moves in this direction

    1. if `lag == 0` new direction is calculated dependent on the social force experienced and an added random error of +-10 degrees; Distance traveled depends on the `speed_factor` (see: `double get_speed`)Movement can occur across window boundaries.
    2. else, new direction is previous direction with added random error but without social force dependency

  #### `double social():`

  - social response is calculated base on the absence or presence of other individuals `in_zone(15)`,`in_zone(100)` and `in_zone(200)` and new direction calculated using `averageturn()`

    1. if the other individuals are within front field of view (330 degrees / 2) and comfort zone (15) focal individual can express:

      1. avoidance (left/right or directly in front) [weight function `(1 / ( 0.1 * pow((dist),2) + 2 ))`]

    2. if the other individuals are within front field of view (330 degrees / 2) and alignment distance (100) focal individual can express:

      1. alignment to closest neighbor [weight function `((1 / ( 0.004 * pow((dist-15),2) + 2 ))`]

    3. if the other individuals are within front field of view (330 degrees / 2) and distance (200) focal individual can express:

      1. attraction to closest neighbor [weight function `(1 / ( 0.004 * pow(( dist-200),2) + 2 ))`]
      2. alignment to closest neighbor [weight function `((1 / ( 0.004 * pow((dist-15),2) + 2 ))`]

  #### `double get_speed():`

  - the `speed_factor` [default=1] is calculated using the individuals in the front visual field (the number of individuals in front (-90,+90 degrees from head tail axis)) and the response radius of 100: speed_factor = 1 + ind_count divided by N/2 (N: Total Population). The more individuals are perceived ahead, the faster the individual will travel.

  <img src="https://github.com/fritzfrancisco/fish_abm/blob/master/Social%20heat%20map.png" width="300">

  Image 2: Social force heat-map

  All visualization was done using OpenCV2 v. 2.4.13 in C++.
