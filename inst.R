# Code to perform instrumental variables analysis comparing surgery with ``PET'' (no surgery)
# Aiming to replicate work by vd Water et al (2012) for the Netherlands

# Number of strategies; 1. unrestricted - just surgery vs no surgery
#                       2. unrestricted, adjust for stage differences
#                       3. unrestricted, fully adjusted
#                       4. restricted based on complete data - 17129 patients used by Jenna
#                       5. restricted complete, fully adjusted
#                       6. restricted after multiple imputation, fully adjusted

# This file is for unrestricted analyses (1-3)

# Step 0. Packages and data

require(survival)
require(dplyr)

regall <- read.csv("regall.csv")

# Step 1. tabulate trusts

# not run
# sort(table(regall$Trust),decreasing = TRUE) # This shows that 33 trusts treated 200 patients or more

gt_200 <- names(which(table(regall$Trust)>200)) # extract these trusts
reg_inst <- regall %>%
              mutate(trust_gt200 = Trust %in% gt_200) %>%
              filter(trust_gt200) %>%
              mutate(Trust = droplevels(Trust)) %>%
              mutate(Trust = factor(Trust, levels = names(sort(table(Trust),decreasing = T))))

# check
# table(reg_inst$Trust)

# Step 2. Create survival objects

reg_inst$status <- 0
reg_inst$status[which(reg_inst$CauseOfDeath == "Breast Cancer")] <- 1
reg_inst$status[which(reg_inst$CauseOfDeath != "Breast Cancer")] <- 2

ri_SurvOS <- with(reg_inst, Surv(Survival, VitalStatus))
ri_SurvBC <- with(reg_inst, Surv(Survival, status == 1))
ri_SurvOC <- with(reg_inst, Surv(Survival, status == 2))

# Step 3. Coxph with Trust

ri_os_coxph1 <- coxph(ri_SurvOS ~  Trust, data = reg_inst)
