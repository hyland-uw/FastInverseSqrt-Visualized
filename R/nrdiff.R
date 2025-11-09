
x_in <- frsr_sample(n = 30000, x_min = 2^-7, x_max = 2.0)$input
three_iters <- filter(frsr(x = x_in, magic = 0x5F400000,
                           tol = 2^-100, detail = TRUE, NRmax = 10), iters == 3)

three_iters <- three_iters %>%
  mutate(one = after_one,
         two = final + diff,
         three = final,
         reference = 1 / sqrt(input)) %>%
  select(input, reference, initial, one, two, three)

three_iters <- three_iters %>%
  mutate(
    init_one = (initial - one) / one,
    one_two = (one - two) / two,
    two_three = (two - three) / three,
    initial = (initial - reference) / reference,
    one = (one - reference)/reference,
    two = (two - reference)/reference,
    three = (three - reference)/reference
  ) %>%
  select(input, initial, init_one, one, one_two, two, two_three, three)



ggplot(three_iters) + geom_point(aes(x = (initial - reference) / reference, y =  one_two / two_three , color = factor(input)), size = 0.4) + coord_polar(theta = "x") + theme_void() + guides(color = "none")



three_long <- three_iters %>%
  pivot_longer(
    cols = c(initial, init_one, one, one_two, two, two_three, three),
    names_to = "variable",
    values_to = "value"
  )
