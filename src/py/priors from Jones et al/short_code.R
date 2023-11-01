




#import original model fit objects
ogawa_rds <- read_rds("source/final_code/luminex_recommendation/generated_rds/decay_model_results/fits-compare/compare-fit-Vibriocidal_Ogawa.rds")
inaba_rds <- read_rds("source/final_code/luminex_recommendation/generated_rds/decay_model_results/fits-compare/compare-fit-Vibriocidal_Inaba.rds")

# request from joseph
# Baseline titers, uninfected (mu_omega)
# Titer increase from baseline upon infection (mu_lambda)
# Delay from reporting to maximum titer rise (mu_logD)
# Half life of titer decay (halflife)


#make draws from posterior
ogawa_draw1000 <- ogawa_rds$fits$exponential_5$fit %>% 
        tidybayes::spread_draws(mu_omega,mu_logD,mu_lambda,halflife,
                                                             ndraws = 1000
                                                             )
inaba_draw1000 <- inaba_rds$fits$exponential_5$fit %>% 
        tidybayes::spread_draws(mu_omega,mu_logD,mu_lambda,halflife,
                                ndraws = 1000
        )

#save objects
write_rds(ogawa_draw1000,"pos_draws_ogawa.rds")
write_rds(inaba_draw1000,"pos_draws_inaba.rds")
write_rds(ogawa_rds$fits$exponential_5,
          "model_obj_ogawa.rds"
          )
write_rds(inaba_rds$fits$exponential_5,
          "model_obj_inaba.rds"
)











