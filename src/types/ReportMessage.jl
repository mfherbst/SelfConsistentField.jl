mutable struct ReportMessage
    msg::String
    data::Any
end
ReportMessage(msg::String) = ReportMessage(msg, nothing)
