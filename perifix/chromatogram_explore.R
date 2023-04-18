raw_chr <- msnexp_filled %>% filterMz(mz=pmppm(guanine_mz, 10)) %>% filterFile(1)
suppressWarnings(plot(raw_chr))


msdata <- msnexp_filled %>%
  fileNames() %>%
  str_subset("Poo") %>%
  grabMSdata()
gp <- msdata$MS1[mz%between%pmppm(guanine_mz)] %>%
  ggplot() +
  geom_line(aes(x=rt, y=int, group=filename))

gsum <- guanine_data %>%
  group_by(feat_id) %>%
  summarize(mz=mean(feat_mzmed), rt=mean(feat_rtmed)/60)
gp + 
  geom_vline(xintercept=gsum$rt, color="red") + 
  ylim(c(0, 1e5)) 
