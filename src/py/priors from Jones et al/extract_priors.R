#library("stan")

model_obj_ogawa <- readRDS("~/Documents/Research/cholera-serology-haiti/analysis_2022/priors from Jones et al/model_obj_ogawa.rds")
ogawa_sum <- model_obj_ogawa$summary
model_obj_inaba <- readRDS("~/Documents/Research/cholera-serology-haiti/analysis_2022/priors from Jones et al/model_obj_inaba.rds")
inaba_sum <- model_obj_inaba$summary

# interesting parameters
print(inaba_sum %>% filter(!str_detect(parameter, '\\[')) %>% select(parameter))

print(paste("(ogawa, inaba) mu_lambda: ", 
            ogawa_sum %>% filter(parameter=="mu_lambda") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="mu_lambda") %>% select("mean")))
print(paste("(ogawa, inaba) mu_omega: ", 
            ogawa_sum %>% filter(parameter=="mu_omega") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="mu_omega") %>% select("mean")))

print(paste("(ogawa, inaba) halflife: ", 
            ogawa_sum %>% filter(parameter=="halflife") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="halflife") %>% select("mean")))

print(paste("(ogawa, inaba) sigma[1,1]: ", 
            ogawa_sum %>% filter(parameter=="params_sigma[1,1]") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="params_sigma[1,1]") %>% select("mean")))
print(paste("(ogawa, inaba) sigma[2,2]: ", 
            ogawa_sum %>% filter(parameter=="params_sigma[2,2]") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="params_sigma[2,2]") %>% select("mean")))
print(paste("(ogawa, inaba) sigma[1,2]: ", 
            ogawa_sum %>% filter(parameter=="params_sigma[1,2]") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="params_sigma[1,2]") %>% select("mean")))
print(paste("(ogawa, inaba) sigma[2,1]: ", 
            ogawa_sum %>% filter(parameter=="params_sigma[2,1]") %>% select("mean"), 
            inaba_sum %>% filter(parameter=="params_sigma[2,1]") %>% select("mean")))
