
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(stringr)
library(DT)



#header_title_javascript = '$(document).ready(function() { $("header").find("nav").append(\'<div id="pageHeader" class="myClass"></div>\');} )'
header_title_javascript = '$(document).ready(function() { $("header").find("nav").append(\'<div id="pageHeader" style="font-size:12px;" class="myClass"></div>\');} )'
collapsa_box_jscode <- "shinyjs.collapse = function(boxid) {$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();}"

#version = "app.v1.4.3"
version = "app.v1.5"
#print(version)
source(file = file.path(version, "ui_app.R"),
         local = TRUE,
         encoding = "UTF-8")




