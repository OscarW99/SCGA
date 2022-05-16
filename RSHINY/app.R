library(shiny)

ui <- fluidPage(
    fluidRow(
        column(
            4,
            textInput(inputId = "path_to_count_data", label = "Provide full path to paient count data")
        )
    )
)

server <- function(input, output) {
}

shinyApp(ui = ui, server = server)