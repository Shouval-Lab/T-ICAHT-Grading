# T-ICAHT coding

## Data preparation

To assign T-ICAHT grades based on longitudinal platelet lab data, please
prepare the following:

1.  A data frame called `df_meta` , which consists of the following
    columns:

<!-- -->

1.  Patient ID (`record_id`)
2.  Time-to-PFS in days after infusion (`tt_pfs_d`)

<!-- -->

2.  A data frame called `df_plt`, which consists of the following
    columns:

<!-- -->

1.  Patient ID (`record_id`)
2.  PLT values (`plt`)
3.  Days relative to CAR-T infusion (`day_rel_car`)

Note that all PLT values should be standardized to the same unit of
(k/mcl).

## Dynamic Early T-ICAHT grading

First, our goal is to create a long format table of patient PLT values
for each of day 0 to day 37. We choose to extend the range from day 30
to day 37 to give more robust grades between day 0 and day 30.

The following code create the long format table, where the PLT
observation is censored at either day 37 or the time of progression.
Missing PLT values will be `NA`s.

``` r
ddist <- sapply(0:37,function(x) as.numeric(df_plt$day_rel_car == x))
colnames(ddist) <- paste("Day",0:37)
ddist <- data.frame(ddist,check.names = FALSE)
ddist$record_id <- df_plt$record_id
ddist$value <- df_plt$plt

ddistlong <- ddist %>%
  pivot_longer(cols=`Day 0`:`Day 37`,
               names_to = "name",
               values_to = "keep") |>
  mutate(day = parse_number(name)) |>
  group_by(record_id,day) |>
  mutate(sumkeep = sum(keep)) |>
  filter((sumkeep == 1 & keep == 1) | (sumkeep == 0 & row_number() == 1)) |> 
  arrange(record_id,day)

ddistlong$value[ddistlong$keep==0] = NA

ddistlong <- ddistlong |> 
  left_join(df_meta |> select(record_id, tt_pfs_d)) |> 
  filter(day <= tt_pfs_d)
```

Next, we perform a simple imputation for the missing PLT values. The
method is simplistic, such that any missing values between two PLT
observations will be assigned as the average of the two observed PLT
values. For missing values before the first (or after the last) observed
PLT value, they will be imputed as the first (or last) value.

``` r
ddistlong_up <- ddistlong |>
  ungroup() |>
  group_by(record_id) |>
  fill(value,.direction="up")

ddistlong_down <- ddistlong |>
  ungroup() |>
  group_by(record_id) |>
  fill(value,.direction="down")

ddistlong_up$value[is.na(ddistlong_up$value)] = ddistlong_down$value[is.na(ddistlong_up$value)]
ddistlong_down$value[is.na(ddistlong_down$value)] = ddistlong_up$value[is.na(ddistlong_down$value)]
ddistlong$value <- (ddistlong_up$value + ddistlong_down$value)/2
ddistlong <- ddistlong[!is.na(ddistlong$value),]
```

Next, we calculate streaks of days when patients with PLT \<= 50.
Specifically, for each row in our long format we annotate the PLT value
as either \<= 50 (`count=1`) or \> 50 (`count=0`). Then we count the
number of consecutive days (`streak1`) when the patient had PLT \<= 50,
as well as the total number of days in each streak (`length_streak`).
Finally, we define the interval between two streaks of PLT \<= 50 as
`inter`.

``` r
df_ticaht_50 <- ddistlong |>
  ungroup() |>
  select(-keep,-sumkeep) |>
  filter(day >= 0 & day <= 37) |>
  group_by(record_id) |>
  mutate(count = as.numeric(value <= 50),
         x = cumsum(c(0,diff(count)) != 0)) |>
  group_by(record_id,x) |>
  mutate(streak1 = ifelse(count==1,row_number(),0),
         length_streak = n()) |>
  group_by(record_id) |>
  mutate(sum_streak1 = cumsum(count==1),
         inter = streak1==0 & sum_streak1 > 0 & sum_streak1 < max(sum_streak1))
```

Up to here, the dynamic T-ICAHT grading is almost done. However, we
would like to add an additional criterion, which only consider the
interval between two streaks an interval if the length of the interval
was greater than 6 days. In other words, for those intervals \<= 6 days,
they were considered too short to separate the two streaks, such that we
consider the two streaks plus the interval between the two streaks a
single long streak.

``` r
df_ticaht$count[df_ticaht$inter & df_ticaht$length_streak <= 6] <- 1
```

With the dynamic status of PLT \<= 50 defined above, it is easy to
define early T-ICAHT grades based on the threshold of 50, where the
variable names is `severe`.

``` r
df_ticaht_below50 <- df_ticaht |>
  select(record_id,day,value,count,center) |>
  group_by(record_id) |>
  mutate(x = cumsum(c(0,diff(count)) != 0)) |>
  group_by(record_id,x) |>
  mutate(streak0 = ifelse(count==0,row_number(),0),
         streak1 = ifelse(count==1,row_number(),0),
         length_streak = n()) |>
  group_by(record_id) |>
  mutate(severe = case_when(streak1 == 0 ~ 0,
                            streak1 < 7 ~ 1,
                            streak1 >= 7 ~ 2))
```

Similarly, we could repeat the above for the other threshold of 20 to
get the early T-ICAHT grade based on the threshold of 20 (the variable
name is `profound`).

``` r
df_ticaht_20 <- ddistlong |>
  ungroup() |>
  select(-keep,-sumkeep) |>
  filter(day >= 0 & day <= 37) |>
  group_by(record_id) |>
  mutate(count = as.numeric(value <= 20),
         x = cumsum(c(0,diff(count)) != 0)) |>
  group_by(record_id,x) |>
  mutate(streak1 = ifelse(count==1,row_number(),0),
         length_streak = n()) |>
  group_by(record_id) |>
  mutate(sum_streak1 = cumsum(count==1),
         inter = streak1==0 & sum_streak1 > 0 & sum_streak1 < max(sum_streak1))

df_ticaht$count[df_ticaht$inter & df_ticaht$length_streak <= 6] <- 1

df_ticaht_below20 <- df_ticaht |>
  select(record_id,day,value,count,center) |>
  group_by(record_id) |>
  mutate(x = cumsum(c(0,diff(count)) != 0)) |>
  group_by(record_id,x) |>
  mutate(streak0 = ifelse(count==0,row_number(),0),
         streak1 = ifelse(count==1,row_number(),0),
         length_streak = n()) |>
  group_by(record_id) |>
  mutate(profound = case_when(streak1 == 0 ~ 0,
                              streak1 < 14 ~ 3,
                              streak1 >= 14 ~ 4))
```

Finally, we could merge the grades based on the 50 criterion and 20
criterion, which forms the final early T-ICAHT dynamic grades
(`ticaht`).

``` r
df_ticaht_combined <- df_ticaht_below50 |>
  left_join(df_ticaht_below20 |> select(record_id, day, profound)) |>
  ungroup() |>
  mutate(ticaht = ifelse(severe > profound, severe, profound))
```

## Late T-ICAHT grading

As specified in the late T-ICAHT definition, we still need the duration
of PLT \<= 20 as a criterion to assign grade 4. Therefore, we applied a
similar coding strategy as we did for early T-ICAHT to create the long
format PLT table between day 31 and day 100. Due to the limited lab
samples beyond day 30 after infusion, we only consider patients with at
least 2 PLT values before the time of progression.

``` r
ddist <- sapply(31:100,function(x) as.numeric(df_plt$day_rel_car == x))
colnames(ddist) <- paste("Day",31:100)
ddist <- data.frame(ddist,check.names = FALSE)
ddist$record_id <- df_plt$record_id
ddist$value <- df_plt$plt

ddistlong <- ddist %>%
  pivot_longer(cols=`Day 31`:`Day 100`,
               names_to = "name",
               values_to = "keep") |>
  mutate(day = parse_number(name)) |>
  group_by(record_id,day) |>
  mutate(sumkeep = sum(keep)) |>
  filter((sumkeep == 1 & keep == 1) | (sumkeep == 0 & row_number() == 1)) |> 
  arrange(record_id,day)

ddistlong$value[ddistlong$keep==0] = NA

ddistlong <- ddistlong |> 
  group_by(record_id) |> 
  mutate(num_obs = sum(!is.na(value))) |> 
  filter(num_obs >= 2)|> 
  left_join(df_meta |> select(record_id, tt_pfs_d)) |> 
  filter(day <= tt_pfs_d) 
```

Based on the late T-ICAHT definition, we assign a ICAHT grade 1, 2, or 3
if there were at least two observations of PLT values below the
corresponding thresholds.

``` r
df_late_lowgrade <- ddistlong[!is.na(ddistlong$value),] |> 
  mutate(grade = case_when(sum(value <= 20) >= 2 ~ 3,
                           sum(value <= 50) >= 2  ~ 2,
                           sum(value <= 100) >= 2 ~ 1,
                           sum(value <= 100) < 2 ~ 0)) |> 
  filter(row_number()==1) |> 
  select(record_id,grade)
```

For grade 4 assignment, we repeat the imputation strategy used for early
T-ICAHT. The difference is that we do not impute the missing values
before the first or after the last PLT observation.

``` r
ddistlong_up <- ddistlong |>
  ungroup() |>
  group_by(record_id) |> 
  fill(value,.direction="up")

ddistlong_down <- ddistlong |>
  ungroup() |>
  group_by(record_id) |>
  fill(value,.direction="down")

ddistlong$value <- (ddistlong_up$value + ddistlong_down$value)/2
ddistlong <- ddistlong[!is.na(ddistlong$value),]

df_late_ticaht <- ddistlong |>
  ungroup() |>
  select(-keep,-sumkeep) |>
  filter(day >= 31 & day <= 100) |>
  group_by(record_id) |>
  mutate(count = as.numeric(value <= 20),
         x = cumsum(c(0,diff(count)) != 0)) |>
  group_by(record_id,x) |>
  mutate(streak1 = ifelse(count==1,row_number(),0),
         length_streak = n()) |>
  group_by(record_id) |>
  mutate(sum_streak1 = cumsum(count==1),
         inter = streak1==0 & sum_streak1 > 0 & sum_streak1 < max(sum_streak1))

df_late_ticaht$count[df_late_ticaht$inter & df_late_ticaht$length_streak <= 6] <- 1

df_late_ticaht_below20 <- df_late_ticaht |>
  select(record_id,day,value,count,center) |>
  group_by(record_id) |>
  mutate(x = cumsum(c(0,diff(count)) != 0)) |>
  group_by(record_id,x) |>
  mutate(streak0 = ifelse(count==0,row_number(),0),
         streak1 = ifelse(count==1,row_number(),0),
         length_streak = n()) |>
  group_by(record_id) |>
  mutate(grade4 = case_when(streak1 == 0 ~ 0,
                            streak1 <= 13 ~ 0,
                            streak1 >= 14 ~ 1)
  )

df_late_icaht <- df_late_ticaht_below20 |>
  group_by(record_id) |> 
  mutate(max_late_ticaht = max(grade4)) |> 
  filter(row_number() == 1) |> 
  ungroup() |> 
  select(record_id,max_late_ticaht) |> 
  left_join(df_late_lowgrade) |> 
  mutate(max_late_ticaht = case_when(max_late_ticaht == 1 ~ 4,
                                     max_late_ticaht == 0 ~ grade))
```

## Summarizing the early and late T-ICAHT grades

Finally, we summarize the early and late T-ICAHT grades by merging the
data into a same data frame.

``` r
df_ticaht_max <- df_ticaht_combined |> 
  group_by(record_id) |> 
  mutate(max_early_ticaht = max(ticaht)) |> 
  filter(row_number() == 1) |> 
  ungroup() |> 
  select(-ticaht) |> 
  left_join(df_late_ticaht |> select(record_id, max_late_ticaht)) |> 
  select(record_id, max_early_ticaht, max_late_ticaht)
```
