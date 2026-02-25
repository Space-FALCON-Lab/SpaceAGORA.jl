output "artifacts_bucket" {
  value = aws_s3_bucket.artifacts.bucket
}

output "run_table" {
  value = aws_dynamodb_table.runs.name
}

output "job_queue_url" {
  value = aws_sqs_queue.jobs.url
}

output "api_endpoint" {
  value = aws_apigatewayv2_stage.default.invoke_url
}
