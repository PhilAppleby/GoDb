package main

import (
	"bufio"
	"ehrdb"
	"fmt"
	"html/template"
	"log"
	"net/http"
	"os"
	"strings"
)

// NONE Substitute for string value "None"
// Found throughout the package
const NONE = "None"

// Convenience function for printing to stdout
func p(a ...interface{}) {
	fmt.Println(a...)
}

// Convenience function to redirect to the error message page
func errorMessage(writer http.ResponseWriter, request *http.Request, msg string) {
	url := []string{"/err?msg=", msg}
	http.Redirect(writer, request, strings.Join(url, ""), 302)
}

// pass in a list of file names, and get a template
func parseTemplateFiles(filenames ...string) (t *template.Template) {
	var files []string
	t = template.New("layout")
	for _, file := range filenames {
		files = append(files, fmt.Sprintf("templates/%s.html", file))
	}
	t = template.Must(t.ParseFiles(files...))
	return
}

func generateHTML(writer http.ResponseWriter, data interface{}, filenames ...string) {
	var files []string
	for _, file := range filenames {
		files = append(files, fmt.Sprintf("templates/%s.html", file))
	}

	templates := template.Must(template.ParseFiles(files...))
	templates.ExecuteTemplate(writer, "layout", data)
}

// for logging
func info(args ...interface{}) {
	logger.SetPrefix("INFO ")
	logger.Println(args...)
}

func danger(args ...interface{}) {
	logger.SetPrefix("ERROR ")
	logger.Println(args...)
}

func warning(args ...interface{}) {
	logger.SetPrefix("WARNING ")
	logger.Println(args...)
}

// version
func version() string {
	return "0.5"
}

func writeBufferedFile(filePath string, records []string) {
	f, err := os.Create(filePath)
	check(err, "Utils, writeBufferedFile (1)")
	defer f.Close()
	w := bufio.NewWriter(f)

	for _, line := range records {
		_, err := w.WriteString(line + "\n")
		check(err, "Utils, writeBufferedFile (2)")
	}
	w.Flush()
}

func convertStringMapToCSV(input map[string]string) []string {
	var lines []string
	lines = make([]string, 0, len(input))
	for k, v := range input {
		lines = append(lines, k+","+v)
	}
	return lines
}

func assocReformat(assocResult string) []string {
	return strings.Split(assocResult, "\n")
}

func getVariantList(urlParams map[string][]string) (string, []string) {
	log.Printf("getVariantList: %v", urlParams)
	varlist := strings.Split(urlParams["variant"][0], ",")
	varlistName := urlParams["varlistname"][0]
	if len(varlist[0]) == 0 {
		if varlistName != NONE {
			varlist, _ = ehrdb.GetVarlistByName(varlistName)
		}
	} else {
		varlistName = NONE
	}
	return varlistName, varlist
}
