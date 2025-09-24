# Carregar pacotes necessários
library(did)
library(dplyr)
library(ggplot2)
library(haven)

# Definir diretório de trabalho
setwd("SUA PASTA AQUI")

# Carregar dados (assumindo que o dataset está em formato .dta)
data <- read_dta("dataset.dta") %>%
  filter(sample1 == 1)

# Criando grupos de tratamento
data <- data %>%
  group_by(state_r) %>%
  mutate(gvar_CS = ifelse(any(L0 == 1), yq[L0 == 1][1], 0)) %>%
  ungroup() %>%
  mutate(gvar_CS = ifelse(is.na(gvar_CS), 0, gvar_CS))

# Lista de variáveis de interesse
dvlist <- c("lnspend_pc", "lnpriceperpresc", "lnpresc_pc", "lnquantity_pcp", 
            "genpen_presc", "geneff_presc", "genacc_presc", "genacc_sim", 
            "s_MCO", "s_offset")



# Função para executar CSDiD e criar gráficos
run_csdid_analysis <- function(depvar, data) {
  
  # CSDiD never treated (controle = nunca tratado)
  cs_nt <- att_gt(yname = depvar,
                  gname = "gvar_CS",
                  idname = "state_r", 
                  tname = "yq",
                  xformla = ~1,
                  data = data,
                  weightsname = "wgt",
                  est_method = "reg",
                  control_group = "nevertreated")
  
  # Event study para never treated
  es_nt <- aggte(cs_nt, type = "dynamic", min_e = -9, max_e = 9, na.rm = TRUE)
  
  # CSDiD not yet treated (controle = ainda não tratado)
  cs_nyt <- att_gt(yname = depvar,
                   gname = "gvar_CS", 
                   idname = "state_r",
                   tname = "yq",
                   xformla = ~1,
                   data = data,
                   weightsname = "wgt",
                   est_method = "reg",
                   control_group = "notyettreated")
  
  # Event study para not yet treated
  es_nyt <- aggte(cs_nyt, type = "dynamic", min_e = -9, max_e = 9, na.rm = TRUE)
  
  # Criar gráfico combinado
  p1 <- ggdid(es_nt, ylim = c(-0.5, 0.5)) + 
    ggtitle("CSDiD Never Treated") +
    theme_minimal()
  
  p2 <- ggdid(es_nyt, ylim = c(-0.5, 0.5)) + 
    ggtitle("CSDiD Not Yet Treated") +
    theme_minimal()
  
  # Combinar dados para gráfico único
  df_nt <- data.frame(
    event_time = es_nt$egt,
    att = es_nt$att.egt,
    se = es_nt$se.egt,
    type = "Never Treated"
  )
  
  df_nyt <- data.frame(
    event_time = es_nyt$egt,
    att = es_nyt$att.egt, 
    se = es_nyt$se.egt,
    type = "Not Yet Treated"
  )
  
  df_combined <- rbind(df_nt, df_nyt) %>%
    mutate(
      ci_lower = att - 1.96 * se,
      ci_upper = att + 1.96 * se
    )
  
  # Gráfico combinado estilo event study
  p_combined <- ggplot(df_combined, aes(x = event_time, y = att, color = type)) +
    geom_point(aes(shape = type), size = 2, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, position = position_dodge(width = 0.3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray") +
    scale_x_continuous(breaks = seq(-9, 9, 1)) +
    scale_color_manual(values = c("Never Treated" = "green", "Not Yet Treated" = "orange")) +
    scale_shape_manual(values = c("Never Treated" = 4, "Not Yet Treated" = 17)) +
    labs(
      title = paste("Event Study -", depvar),
      x = "Event Time",
      y = "ATT",
      color = "Control Group",
      shape = "Control Group"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  print(p_combined)
  
  # Salvar gráfico
  ggsave(paste0("replication_results/id_51 - R - ", depvar, ".pdf"), 
         plot = p_combined, width = 10, height = 6)
  
  # Retornar resultados
  return(list(
    never_treated = list(cs = cs_nt, es = es_nt),
    not_yet_treated = list(cs = cs_nyt, es = es_nyt),
    plot = p_combined
  ))
}

# Criar diretório para resultados se não existir
dir.create("replication_results", showWarnings = FALSE)

# Loop para estimar o DiD para todas as variáveis
results_list <- list()

for(depvar in dvlist) {
  cat("Processando variável:", depvar, "\n")
  
  # Executar análise
  results_list[[depvar]] <- run_csdid_analysis(depvar, data)
  
  # Imprimir resumo dos resultados
  cat("\n=== Resultados para", depvar, "===\n")
  cat("Never Treated:\n")
  print(summary(results_list[[depvar]]$never_treated$es))
  cat("\nNot Yet Treated:\n") 
  print(summary(results_list[[depvar]]$not_yet_treated$es))
  cat("\n", paste(rep("=", 50), collapse = ""), "\n\n")
}

